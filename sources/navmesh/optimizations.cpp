#include <atomic>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "navmesh.h"

#include <cage-core/concurrent.h>
#include <cage-core/tasks.h>

namespace std
{
	template<>
	struct hash<std::pair<cage::uint32, cage::uint32>>
	{
		size_t operator()(const std::pair<cage::uint32, cage::uint32> &x) const { return hash<cage::uint32>()(x.first) ^ hash<cage::uint32>()(x.first ^ x.second); }
	};
}

namespace navoptim
{
	constexpr Real shortEdge = 0.55f;
	constexpr Real longEdge = 1.3f;
	constexpr Real longEdgeThreshold = 0.6f;
	static_assert(longEdge > longEdgeThreshold * 2);
	static_assert(longEdgeThreshold > shortEdge);

	std::vector<uint32> sharedNeighbors(const Graph &graph, uint32 a, uint32 b)
	{
		const auto &ma = graph.neighbors.at(a);
		const auto &mb = graph.neighbors.at(b);
		std::vector<uint32> r;
		r.resize(min(ma.size(), mb.size()));
		r.resize(std::set_intersection(ma.begin(), ma.end(), mb.begin(), mb.end(), r.begin()) - r.begin());
		return r;
	}

	FlatSet<uint32> connectedNodes(const Graph &graph, const FlatSet<uint32> &start)
	{
		FlatSet<uint32> res;
		res.reserve(start.size() * 2 + 10);
		for (uint32 a : start)
			for (uint32 b : graph.neighbors[a])
				res.insert(b);
		return res;
	}

	uint32 totalEdgesCount(const Graph &graph)
	{
		uint32 sum = 0;
		for (const auto &it : graph.neighbors)
			sum += numeric_cast<uint32>(it.size());
		return sum / 2;
	}

	void removeEmptyNodes(Graph &graph)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "removing empty nodes");
		const uint32 initialNodesCount = numeric_cast<uint32>(graph.nodes.size());

		std::vector<bool> deleting;
		deleting.reserve(initialNodesCount);
		for (const auto &it : graph.neighbors)
			deleting.push_back(it.empty());

		{
			std::vector<uint32> remap;
			remap.reserve(initialNodesCount);
			for (bool d : deleting)
				remap.push_back(d);
			for (uint32 i = 1; i < initialNodesCount; i++)
				remap[i] += remap[i - 1];
			for (uint32 i = 0; i < initialNodesCount; i++)
			{
				for (uint32 &b : graph.neighbors[i].unsafeData())
					b -= remap[b];
			}
		}

		{
			auto &ns = graph.nodes;
			auto it = deleting.begin();
			std::erase_if(ns, [&](const Node &) { return *it++; });
		}
		{
			auto &ns = graph.neighbors;
			auto it = deleting.begin();
			std::erase_if(ns, [&](const FlatSet<uint32> &) { return *it++; });
		}

		graphValidationDebugOnly(graph);
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "nodes count: " + graph.nodes.size() + " (was " + initialNodesCount + ")");
	}

	void markBorderNodes(Graph &graph)
	{
		// this MUST be done before the additional counter-diagonal links are added
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "marking border vertices");

		const auto &border = [&](uint32 a)
		{
			if (graph.neighbors[a].size() < 3)
				return true;
			for (uint32 b : graph.neighbors[a])
			{
				if (sharedNeighbors(graph, a, b).size() < 2)
					return true;
			}
			return false;
		};
		uint32 index = 0;
		for (auto &it : enumerate(graph.nodes))
		{
			it->border = border(numeric_cast<uint32>(it.index));
			if (it->border)
				index++;
		}

		graphValidationDebugOnly(graph);
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "borders count: " + index);
	}

	void optimizeNodePositions(Graph &graph, Real factor)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "optimizing node positions");

		struct Runner : Immovable
		{
			Graph &graph;
			std::vector<Vec3> forces;
			std::vector<FlatSet<uint32>> nearby;
			const uint32 count = numeric_cast<uint32>(graph.nodes.size());
			Real factor;

			Runner(Graph &graph, Real factor) : graph(graph), factor(factor) {}

			static Real springForce(Real x)
			{
				const auto &a = [](Real x) { return sin(Rads(-Real::Pi() * (x * 1.9 + 0.5))); };
				const auto &b = [](Real x) { return 1 / (x + 0.1); };
				return a(x) * b(x);
			};

			void operator()(uint32 invocation)
			{
				const auto range = tasksSplit(invocation, processorsCount(), count);
				for (uint32 a = range.first; a < range.second; a++)
				{
					if (graph.nodes[a].border)
						continue; // border nodes do not move

					const Vec3 &ap = graph.nodes[a].position;
					Vec3 &force = forces[a];
					force = Vec3();

					// spring forces
					for (uint32 b : nearby[a])
					{
						const Vec3 &bp = graph.nodes[b].position;
						const Vec3 d = ap - bp;
						const Real l = length(d);
						const Real s = springForce(l);
						force += normalize(d) * s;
					}

					// project the force into the plane of the node
					const Vec3 &n = graph.nodes[a].normal;
					force -= n * dot(force, n);
				}
			}

			void run()
			{
				forces.resize(count);
				nearby.resize(count);
				for (const auto &it : enumerate(graph.nodes))
				{
					FlatSet<uint32> ngs = connectedNodes(graph, { numeric_cast<uint32>(it.index) });
					ngs = connectedNodes(graph, ngs);
					ngs = connectedNodes(graph, ngs);
					ngs.erase(numeric_cast<uint32>(it.index));
					nearby[it.index] = ngs;
				}

				for (uint32 iteration = 0; iteration < 13; iteration++)
				{
					CAGE_LOG_DEBUG(SeverityEnum::Info, "libnavmesh", Stringizer() + "node positions, iteration: " + iteration);

					tasksRunBlocking<Runner>("optimizing node positions", *this, processorsCount());

					// apply the forces
					for (uint32 a = 0; a < count; a++)
						graph.nodes[a].position += forces[a] * 0.015 * factor;
				}
			}
		};

		{
			Runner runner(graph, factor);
			runner.run();
		}

		graphValidationDebugOnly(graph);
	}

	void snapNodesToCollider(Graph &graph, const Collider *collider, Real tileSize)
	{
		if (!collider)
			return;

		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "snapping nodes to collider");

		struct Runner
		{
			Graph &graph;
			const Collider *collider;
			Real tileSize;
			std::atomic<uint32> warnings = 0;

			Runner(Graph &graph, const Collider *collider, Real tileSize) : graph(graph), collider(collider), tileSize(tileSize) {}

			void operator()(uint32 invocation)
			{
				uint32 warn = 0;
				const auto range = tasksSplit(invocation, processorsCount(), numeric_cast<uint32>(graph.nodes.size()));
				for (uint32 index = range.first; index < range.second; index++)
				{
					const Vec3 a = graph.nodes[index].position;
					const Vec3 n = graph.nodes[index].normal;
					const Line l = makeSegment(a - n, a + n);
					const Vec3 p = intersection(l, collider, Transform({}, {}, 1 / tileSize));
					if (valid(p))
						graph.nodes[index].position = p;
					else if (!graph.nodes[index].border) // border nodes may extend slightly outside the collider
						warn++;
				}
				warnings += warn;
			}

			void run() { tasksRunBlocking<Runner>("snapping nodes to collider", *this, processorsCount()); }
		};

		{
			Runner runner(graph, collider, tileSize);
			runner.run();
			if (runner.warnings > 0)
				CAGE_LOG(SeverityEnum::Warning, "libnavmesh", Stringizer() + (uint32)runner.warnings + " nodes could not snap to collider");
		}
	}

	void addNeighborEdges(Graph &graph)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "adding neighbor edges");
		const uint32 initialEdgesCount = totalEdgesCount(graph);

		std::unordered_set<std::pair<uint32, uint32>> additions;
		for (const auto &it : enumerate(graph.neighbors))
		{
			for (uint32 n1 : *it)
			{
				for (uint32 n2 : graph.neighbors[n1])
				{
					if (it.index >= n2)
						continue;
					if (graph.nodes[it.index].border && graph.nodes[n2].border)
						continue; // path going through "outside" of the map
					additions.emplace(numeric_cast<uint32>(it.index), n2);
				}
			}
		}

		for (const auto &it : additions)
		{
			graph.neighbors[it.first].insert(it.second);
			graph.neighbors[it.second].insert(it.first);
		}

		graphValidationDebugOnly(graph);
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "edges count: " + totalEdgesCount(graph) + " (was " + initialEdgesCount + ")");
	}

	void splitLongEdges(Graph &graph)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "splitting long edges");
		const uint32 initialNodesCount = numeric_cast<uint32>(graph.nodes.size());

		std::unordered_set<std::pair<uint32, uint32>> additions;
		for (uint32 a = 0; a < initialNodesCount; a++)
		{
			const Vec3 &ap = graph.nodes[a].position;
			for (uint32 b : graph.neighbors[a])
			{
				if (a >= b)
					continue;
				const Vec3 &bp = graph.nodes[b].position;
				if (distanceSquared(ap, bp) < sqr(longEdge))
					continue;
				const Vec3 mp = (ap + bp) * 0.5;
				bool ok = true;
				for (uint32 c : graph.neighbors[a])
				{
					if (distanceSquared(mp, graph.nodes[c].position) < sqr(longEdgeThreshold))
					{
						ok = false;
						break;
					}
				}
				if (!ok)
					continue;
				for (uint32 c : graph.neighbors[b])
				{
					if (distanceSquared(mp, graph.nodes[c].position) < sqr(longEdgeThreshold))
					{
						ok = false;
						break;
					}
				}
				if (!ok)
					continue;
				additions.emplace(a, b);
			}
		}

		for (const auto &it : additions)
		{
			CAGE_ASSERT(it.first < it.second);
			const uint32 c = numeric_cast<uint32>(graph.nodes.size());
			{
				graph.nodes.emplace_back();
				auto &v = graph.nodes[c];
				const auto &a = graph.nodes[it.first];
				const auto &b = graph.nodes[it.second];
				v.position = (a.position + b.position) * 0.5;
				v.normal = normalize(a.normal + b.normal);
			}
			{
				graph.neighbors.emplace_back();
				auto &v = graph.neighbors[c];
				const auto &a = graph.neighbors[it.first];
				const auto &b = graph.neighbors[it.second];
				v.insert(a.begin(), a.end());
				v.insert(b.begin(), b.end());
				for (uint32 n : a)
					graph.neighbors[n].insert(c);
				for (uint32 n : b)
					graph.neighbors[n].insert(c);
				v.erase(c);
			}
		}

		for (const auto &it : additions)
		{
			graph.neighbors[it.first].erase(it.second);
			graph.neighbors[it.second].erase(it.first);
		}

		graphValidationDebugOnly(graph);
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "nodes count: " + graph.nodes.size() + " (was " + initialNodesCount + ")");
	}

	void joinCloseNodes(Graph &graph)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "joining close nodes");

		// find which nodes to join and where
		std::unordered_map<uint32, uint32> joins;
		for (const auto &it : enumerate(graph.neighbors))
		{
			const uint32 a = numeric_cast<uint32>(it.index);
			const Vec3 &ap = graph.nodes[a].position;
			for (uint32 b : *it)
			{
				if (a >= b)
					continue;
				const Vec3 &bp = graph.nodes[b].position;
				if (distanceSquared(ap, bp) < sqr(shortEdge))
					joins[b] = a;
			}
		}

		// make sure that, when joining more than two nodes, all join onto one
		for (auto &it : joins)
		{
			while (joins.count(it.second) > 0)
				it.second = joins[it.second];
		}

		// all edges going to or from A are replaced to go to or from B
		const auto &moveEdges = [&](uint32 a, uint32 b)
		{
			CAGE_ASSERT(a != b);
			auto &m = graph.neighbors[a];
			auto &t = graph.neighbors[b];
			for (uint32 n : m)
			{
				auto &k = graph.neighbors[n];
				k.erase(a);
				k.insert(b);
			}
			t.insert(m.begin(), m.end());
			t.erase(b);
			m.clear();
		};

		// reconnect vertices
		for (const auto &it : joins)
			moveEdges(it.first, it.second);

		removeEmptyNodes(graph);
	}

	void removeSuboptimalEdges(Graph &graph)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "removing suboptimal edges");
		const uint32 initialEdgesCount = totalEdgesCount(graph);

		struct Runner
		{
			Graph &graph;
			std::unordered_set<std::pair<uint32, uint32>> deletions;
			Holder<Mutex> mutex = newMutex();

			Runner(Graph &graph) : graph(graph) {}

			void operator()(uint32 invocation)
			{
				std::unordered_set<std::pair<uint32, uint32>> dels;
				dels.reserve(graph.neighbors.size());

				const auto range = tasksSplit(invocation, processorsCount(), numeric_cast<uint32>(graph.nodes.size()));
				for (uint32 index = range.first; index < range.second; index++)
				{
					const Vec3 &a = graph.nodes[index].position;
					for (uint32 n1 : graph.neighbors[index])
					{
						const Vec3 &b = graph.nodes[n1].position;
						const Line rb = makeRay(a, b);
						if (!rb.valid())
							continue;
						for (uint32 n2 : graph.neighbors[index])
						{
							if (n1 >= n2)
								continue;
							const Vec3 &c = graph.nodes[n2].position;
							const Line rc = makeRay(a, c);
							if (!rc.valid())
								continue;

							if (dot(rb.direction, rc.direction) > 0.8191520442)
							{
								const uint32 n = distanceSquared(a, b) > distanceSquared(a, c) ? n1 : n2;
								dels.emplace(min(n, index), max(n, index));
							}
						}
					}
				}

				{
					ScopeLock lock(mutex);
					deletions.merge(dels);
				}

				dels.clear();
			}

			void run()
			{
				deletions.reserve(graph.neighbors.size() * 10);

				tasksRunBlocking<Runner>("removing suboptimal edges", *this, processorsCount());

				for (const auto &it : deletions)
				{
					graph.neighbors[it.first].erase(it.second);
					graph.neighbors[it.second].erase(it.first);
				}
			}
		};

		{
			Runner runner(graph);
			runner.run();
		}

		graphValidationDebugOnly(graph);
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "edges count: " + totalEdgesCount(graph) + " (was " + initialEdgesCount + ")");

		removeEmptyNodes(graph);
	}

	void updateNodeProperties(Graph &graph, const SpatialGraph &original)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "updating node properties");
		for (const auto &it : enumerate(graph.nodes))
		{
			original.spatialQuery->intersection(Sphere(it->position, 3));
			const auto &res = original.spatialQuery->result();
			if (res.size() == 0)
			{
				CAGE_LOG_THROW(Stringizer() + "position: " + it->position);
				CAGE_THROW_ERROR(Exception, "navigation node too far from any original vertex");
			}
			uint32 bestIndex = m;
			Real bestDist = Real::Infinity();
			for (const auto &on : res)
			{
				const Real d = distanceSquared(original.nodes.at(on).position, it->position);
				if (d < bestDist)
				{
					bestDist = d;
					bestIndex = on;
				}
			}
			const auto &orig = original.nodes.at(bestIndex);
			it->normal = orig.normal;
			it->terrain = orig.terrain;
		}
	}
}
