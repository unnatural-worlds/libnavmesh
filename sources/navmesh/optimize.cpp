#include <cage-core/geometry.h>
#include <cage-core/spatialStructure.h>

#include "navmesh.h"

#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/SurfaceRemeshing.h>

#include <algorithm>
#include <utility>
#include <cmath>

namespace cage
{
	struct VecComp
	{
		bool operator () (const vec3 &a, const vec3 &b) const
		{
			return detail::memcmp(&a, &b, sizeof(vec3)) < 0;
		}
	};
}

namespace
{
	constexpr float shortEdge = 0.55;
	constexpr float longEdge = 1.3;
	constexpr float longEdgeThreshold = 0.6;

	NavigationData navigation;
	NavigationData origNeighs;
	Holder<SpatialStructure> origSpatData;
	Holder<SpatialQuery> origSpatQuery;

	uint32 countConnections()
	{
		uint32 nt = 0;
		for (const auto &it : navigation)
			nt += numeric_cast<uint32>(it.second.neighbors.size());
		return nt;
	}

	void printNavigationStatistics()
	{
		real lenSum;
		real lenMin = real::Infinity();
		real lenMax;
		uint32 nghCnt = 0;
		for (const auto &it : navigation)
		{
			for (uint32 nid : it.second.neighbors)
			{
				const auto &n = navigation[nid];
				real d = distance(it.second.position, n.position);
				lenSum += d;
				lenMin = min(lenMin, d);
				lenMax = max(lenMax, d);
				nghCnt++;
			}
		}
		real lenAvg = lenSum / nghCnt;
		real devSum1 = 0;
		real devSum2 = 0;
		for (const auto &it : navigation)
		{
			for (uint32 nid : it.second.neighbors)
			{
				const auto &n = navigation[nid];
				real d = distance(it.second.position, n.position);
				devSum1 += abs(d - lenAvg);
				devSum2 += abs(d - 1);
			}
		}
		real devAvg1 = devSum1 / nghCnt;
		real devAvg2 = devSum2 / nghCnt;
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "average edge length: " + lenAvg + ", deviation: " + devAvg1 + ", fitness: " + devAvg2);
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "minimum: " + lenMin + ", maximum: " + lenMax);
	}

	void prepareOrigData()
	{
		origNeighs = navigation;
		SpatialStructureCreateConfig sdcc;
		sdcc.maxItems = 1000000;
		origSpatData = newSpatialStructure(sdcc);
		origSpatQuery = newSpatialQuery(origSpatData.get());
		for (const auto &it : navigation)
			origSpatData->update(it.first, aabb(it.second.position));
		origSpatData->rebuild();
	}

	void clearOrigData()
	{
		origSpatQuery.clear();
		origSpatData.clear();
		origNeighs.clear();
	}

	std::vector<uint32> commonNeighbors(const NavigationData &map, uint32 a, uint32 b)
	{
		std::vector<uint32> r;
		const auto &ma = map.at(a);
		const auto &mb = map.at(b);
		r.resize(min(ma.neighbors.size(), mb.neighbors.size()));
		r.resize(std::set_intersection(ma.neighbors.begin(), ma.neighbors.end(), mb.neighbors.begin(), mb.neighbors.end(), r.begin()) - r.begin());
		return r;
	}

	void validateNavigationRelease()
	{
		for (const auto &it : navigation)
		{
			if (it.first >= navigation.size())
				CAGE_THROW_CRITICAL(Exception, "validation error: vertex indices must be consecutive");
			if (it.second.neighbors.empty())
				CAGE_THROW_CRITICAL(Exception, "validation error: vertex has no edges");
			for (uint32 n : it.second.neighbors)
			{
				if (navigation.count(n) == 0)
					CAGE_THROW_CRITICAL(Exception, "validation error: invalid edge index");
				if (navigation[n].neighbors.count(it.first) != 1)
					CAGE_THROW_CRITICAL(Exception, "validation error: edge is missing counterpart");
				if (it.first == n)
					CAGE_THROW_CRITICAL(Exception, "validation error: edge to itself");
			}
		}
	}

	void validateNavigationDebug()
	{
#ifdef CAGE_DEBUG
		validateNavigationRelease();
#endif // CAGE_DEBUG
	}

	std::set<uint32> getConnected(const std::set<uint32> &start)
	{
		std::set<uint32> res;
		for (uint32 a : start)
			for (uint32 b : navigation[a].neighbors)
				res.insert(b);
		return res;
	}

	real erf(real x)
	{
		return std::erf(x.value);
	};

	std::vector<vec3> getPositions()
	{
		std::vector<vec3> res;
		res.reserve(navigation.size());
		for (const auto &it : navigation)
			res.push_back(it.second.position);
		return res;
	}

	void joinIdenticalVertices()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "joining identical vertices");

		std::map<vec3, uint32, VecComp> vertices;
		for (const auto &it : navigation)
			vertices[it.second.position];
		{ // assign new indices
			uint32 i = 0;
			for (auto &it : vertices)
				it.second = i++;
		}

		NavigationData nav;
		for (const auto &it : navigation)
		{
			uint32 idx = vertices[it.second.position];
			NavigationPoint &p = nav[idx];
			p.position = it.second.position;
			p.normal += it.second.normal;
			p.terrain = max(p.terrain, it.second.terrain); // prefer the more rare terrain (assuming they are sorted)
			for (auto &n : it.second.neighbors)
				p.neighbors.insert(vertices[navigation[n].position]);
		}
		for (auto &n : nav)
			n.second.normal = normalize(n.second.normal);
		std::swap(navigation, nav);

		validateNavigationDebug();
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "vertices count: " + navigation.size() + " (was " + nav.size() + ")");
	}

	void markBorderVertices()
	{
		// this MUST be done before the additional counter-diagonal links are added
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "marking border vertices");
		const auto &border = [&](uint32 a)
		{
			if (navigation.at(a).neighbors.size() < 3)
				return true;
			for (uint32 o : navigation.at(a).neighbors)
			{
				if (commonNeighbors(navigation, a, o).size() < 2)
					return true;
			}
			return false;
		};
		uint32 cnt = 0;
		for (auto &it : navigation)
		{
			it.second.border = border(it.first);
			if (it.second.border)
				cnt++;
		}
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "borders count: " + cnt);
	}

	void addCounterDiagonals()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "adding counter-diagonal edges");
		uint32 cons = countConnections();

		std::set<std::pair<uint32, uint32>> additions;
		const auto &diagonal = [&](uint32 a, uint32 b)
		{
			auto cs = commonNeighbors(navigation, a, b);
			if (cs.size() == 2)
				additions.emplace(cs[0], cs[1]);
		};

		uint32 vertexCount = numeric_cast<uint32>(navigation.size());
		const std::vector<vec3> positions = getPositions();
		for (uint32 a = 0; a < vertexCount; a++)
		{
			const vec3 ap = positions[a];
			for (uint32 b : navigation[a].neighbors)
			{
				const vec3 bp = positions[b];
				for (uint32 c : commonNeighbors(navigation, a, b))
				{
					const vec3 cp = positions[c];
					real ab = distanceSquared(ap, bp);
					real ac = distanceSquared(ap, cp);
					real bc = distanceSquared(bp, cp);
					if (ab > ac && ab > bc)
						diagonal(a, b);
					else if (ac > ab && ac > bc)
						diagonal(a, c);
					else if (bc > ab && bc > ac)
						diagonal(b, c);
				}
			}
		}

		for (const auto &it : additions)
		{
			navigation[it.first].neighbors.insert(it.second);
			navigation[it.second].neighbors.insert(it.first);
		}

		validateNavigationDebug();
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "edges count: " + countConnections() + " (was " + cons + ")");
	}

	void addNeighborEdges()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "adding neighbor edges");
		uint32 cons = countConnections();

		std::set<std::pair<uint32, uint32>> additions;
		for (auto &it : navigation)
		{
			for (uint32 n1 : it.second.neighbors)
			{
				for (uint32 n2 : navigation[n1].neighbors)
				{
					if (it.first >= n2)
						continue;
					const auto &s = navigation[n2];
					if (it.second.border && s.border)
						continue; // path going through "outside" of the map
					additions.emplace(it.first, n2);
				}
			}
		}

		for (const auto &it : additions)
		{
			navigation[it.first].neighbors.insert(it.second);
			navigation[it.second].neighbors.insert(it.first);
		}

		validateNavigationDebug();
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "edges count: " + countConnections() + " (was " + cons + ")");
	}

	void removeSuboptimalEdges()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "removing suboptimal edges");
		uint32 cons = countConnections();

		const std::vector<vec3> positions = getPositions();
		std::set<std::pair<uint32, uint32>> deletions;
		for (const auto &it : navigation)
		{
			const vec3 a = it.second.position;
			for (uint32 n1 : it.second.neighbors)
			{
				const vec3 b = positions[n1];
				line rb = makeRay(a, b);
				if (!rb.valid())
					continue;
				for (uint32 n2 : it.second.neighbors)
				{
					if (n1 >= n2)
						continue;
					const vec3 c = positions[n2];
					line rc = makeRay(a, c);
					if (!rc.valid())
						continue;
					//if (angle(rb, rc) < degs(35))
					if (dot(rb.direction, rc.direction) > 0.8191520442) // cos(degs(35))
					{
						uint32 n = distanceSquared(a, b) > distanceSquared(a, c) ? n1 : n2;
						deletions.emplace(min(it.first, n), max(it.first, n));
					}
				}
			}
		}

		for (const auto &it : deletions)
		{
			navigation[it.first].neighbors.erase(it.second);
			navigation[it.second].neighbors.erase(it.first);
		}

		validateNavigationDebug();
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "edges count: " + countConnections() + " (was " + cons + ")");
	}

	void splitLongEdges()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "splitting long edges");

		uint32 verticesCount = numeric_cast<uint32>(navigation.size());
		const std::vector<vec3> positions = getPositions();
		std::set<std::pair<uint32, uint32>> additions;
		for (uint32 a = 0; a < verticesCount; a++)
		{
			const vec3 ap = positions[a];
			const std::set<uint32> &aNs = navigation[a].neighbors;
			for (uint32 b : aNs)
			{
				if (a >= b)
					continue;
				const vec3 bp = positions[b];
				if (distanceSquared(ap, bp) < sqr(longEdge))
					continue;
				const vec3 mp = (ap + bp) * 0.5;
				bool ok = true;
				for (uint32 c : aNs)
				{
					if (distanceSquared(mp, positions[c]) < sqr(longEdgeThreshold))
					{
						ok = false;
						break;
					}
				}
				if (!ok)
					continue;
				for (uint32 c : navigation[b].neighbors)
				{
					if (distanceSquared(mp, positions[c]) < sqr(longEdgeThreshold))
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
			CAGE_ASSERT(it.first != it.second);
			uint32 c = numeric_cast<uint32>(navigation.size());
			CAGE_ASSERT(navigation.count(c) == 0);
			auto &v = navigation[c];
			auto &a = navigation[it.first];
			auto &b = navigation[it.second];
			v.position = (a.position + b.position) * 0.5;
			v.normal = normalize(a.normal + b.normal);
			v.neighbors.insert(a.neighbors.begin(), a.neighbors.end());
			v.neighbors.insert(b.neighbors.begin(), b.neighbors.end());
			for (uint32 n : a.neighbors)
				navigation[n].neighbors.insert(c);
			for (uint32 n : b.neighbors)
				navigation[n].neighbors.insert(c);
			v.neighbors.erase(c);
		}

		validateNavigationDebug();
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "vertices count: " + navigation.size() + " (was " + verticesCount + ")");
	}

	void joinCloseVertices()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "joining close vertices");
		uint32 verticesCount = numeric_cast<uint32>(navigation.size());

		// find what vertices to join and where
		std::map<uint32, uint32> joins;
		const std::vector<vec3> positions = getPositions();
		for (const auto &it : navigation)
		{
			uint32 a = it.first;
			const vec3 ap = it.second.position;
			for (uint32 b : it.second.neighbors)
			{
				if (a >= b)
					continue;
				const vec3 bp = positions[b];
				if (distanceSquared(ap, bp) < sqr(shortEdge))
					joins[b] = a;
			}
		}

		// make sure that, when joining more than two vertices, all join onto one
		for (auto &it : joins)
		{
			while (joins.count(it.second) > 0)
				it.second = joins[it.second]; // is this guaranteed to finish??
		}

		// all edges going to or from A are replaced to go to or from B
		const auto &moveEdges = [](uint32 a, uint32 b)
		{
			CAGE_ASSERT(a != b);
			auto &m = navigation[a];
			auto &t = navigation[b];
			for (uint32 n : m.neighbors)
			{
				auto &k = navigation[n];
				k.neighbors.erase(a);
				k.neighbors.insert(b);
			}
			t.neighbors.insert(m.neighbors.begin(), m.neighbors.end());
			t.neighbors.erase(b);
			m.neighbors.clear();
		};

		// reconnect vertices
		for (const auto &it : joins)
			moveEdges(it.first, it.second);

		// move empty vertices to end
		{
			uint32 b = numeric_cast<uint32>(navigation.size() - 1);
			for (const auto &it : joins)
			{
				uint32 a = it.first;
				while (joins.count(b) || navigation[b].neighbors.empty())
					b--;
				if (a >= b)
					continue;
				moveEdges(b, a);
				std::swap(navigation[a], navigation[b]);
				std::swap(navigation[a].neighbors, navigation[b].neighbors);
				b--;
			}
		}

		// actually remove the vertices
		uint32 targetCount = numeric_cast<uint32>(navigation.size() - joins.size());
		while (navigation.size() > targetCount)
		{
			CAGE_ASSERT(navigation.rbegin()->second.neighbors.empty());
			navigation.erase(navigation.rbegin()->first);
		}

		validateNavigationDebug();
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "vertices count: " + navigation.size() + " (was " + verticesCount + ")");
	}

	void removeEmptyVertices()
	{
		uint32 i = 0;
		while (i < navigation.size())
		{
			CAGE_ASSERT(navigation.count(i) == 1);
			if (!navigation[i].neighbors.empty())
			{
				i++;
				continue;
			}
			CAGE_LOG(SeverityEnum::Warning, "libnavmesh", stringizer() + "removing empty vertex: " + i);
			if (navigation[i].border)
				CAGE_LOG(SeverityEnum::Warning, "libnavmesh", "the vertex is at border");
			for (uint32 j = i, e = numeric_cast<uint32>(navigation.size()) - 1; j != e; j++)
				navigation[j] = templates::move(navigation[j + 1]);
			navigation.erase(numeric_cast<uint32>(navigation.size()) - 1);
			for (auto &it : navigation)
			{
				const auto ns = templates::move(it.second.neighbors);
				for (uint32 n : ns)
					it.second.neighbors.insert(n < i ? n : n - 1);
			}
		}
	}

	void optimizeVertexPositions()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "optimizing vertex positions");

		std::vector<vec3> originalPositions, normals, positions, forces;
		std::vector<std::vector<uint32>> nearby, neighbors;
		std::vector<bool> borders;
		const uint32 count = numeric_cast<uint32>(navigation.size());
		originalPositions.reserve(count);
		normals.reserve(count);
		positions.reserve(count);
		nearby.resize(count);
		neighbors.resize(count);
		forces.resize(count);
		borders.reserve(count);
		for (const auto &it : navigation)
		{
			originalPositions.push_back(it.second.position);
			normals.push_back(it.second.normal);
			positions.push_back(it.second.position);
			borders.push_back(it.second.border);
			std::set<uint32> ngs = getConnected({ it.first });
			ngs = getConnected(ngs);
			ngs = getConnected(ngs);
			ngs.erase(it.first);
			nearby[it.first] = std::vector<uint32>(ngs.begin(), ngs.end());
			neighbors[it.first] = std::vector<uint32>(it.second.neighbors.begin(), it.second.neighbors.end());
		}

		const auto &springForce = [](real x) -> real
		{
			const auto &a = [](real x) { return sin(rads(-real::Pi() * (x * 1.9 + 0.5))); };
			const auto &b = [](real x) { return 1 / (x + 0.1); };
			return a(x) * b(x);
		};

		for (uint32 iteration = 0; iteration < 13; iteration++)
		{
			CAGE_LOG_DEBUG(SeverityEnum::Info, "libnavmesh", stringizer() + "vertex positions, iteration: " + iteration);

			for (uint32 a = 0; a < count; a++)
			{
				if (borders[a])
					continue; // border vertices may not move

				vec3 ap = positions[a];
				vec3 &force = forces[a];
				force = vec3();

				// spring forces
				for (uint32 b : nearby[a])
				{
					vec3 bp = positions[b];
					vec3 d = ap - bp;
					real l = length(d);
					real s = springForce(l);
					force += normalize(d) * s;
				}

				// project the force into the plane of the vertex
				force -= normals[a] * dot(force, normals[a]);
			}

			for (uint32 a = 0; a < count; a++)
				positions[a] += forces[a] * 0.015;
		}

		// apply changes
		for (auto &it : navigation)
			it.second.position = positions[it.first];

		validateNavigationDebug();
	}

	void updateVertexPropertiesFromOrig()
	{
		for (auto &it : navigation)
		{
			origSpatQuery->intersection(sphere(it.second.position, 3));
			if (origSpatQuery->result().size() == 0)
				CAGE_THROW_ERROR(Exception, "navigation vertex too far from any original vertex");
			uint32 bestIndex = origSpatQuery->result()[0];
			real bestDist = distanceSquared(origNeighs[origSpatQuery->result()[0]].position, it.second.position);
			for (const auto &on : origSpatQuery->result())
			{
				real d = distanceSquared(origNeighs[on].position, it.second.position);
				if (d < bestDist)
				{
					bestDist = d;
					bestIndex = on;
				}
			}
			const auto &orig = origNeighs[bestIndex];
			it.second.normal = orig.normal;
			it.second.terrain = orig.terrain;
		}
	}

	void pmpRegularization()
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "PMP regularization");
		pmp::SurfaceMesh pm;

		{
			// add vertices
			uint32 idx = 0;
			for (const auto &it : navigation)
			{
				CAGE_ASSERT(it.first == idx++); // ensure that input vertices are continuous
				auto r = pm.add_vertex(pmp::Point(it.second.position[0].value, it.second.position[1].value, it.second.position[2].value));
				CAGE_ASSERT(it.first == r.idx()); // ensure that output vertices are continuous
			}
		}
		{
			// add triangles
			for (const auto &a_ : navigation)
			{
				uint32 a = a_.first;
				for (uint32 b : a_.second.neighbors)
				{
					if (b < a)
						continue;
					const auto &b_ = navigation[b];
					for (uint32 c : commonNeighbors(navigation, a, b))
					{
						if (c < b)
							continue;
						const auto &c_ = navigation[c];
						vec3 tn = normalize(a_.second.normal + b_.normal + c_.normal);
						vec3 mn = triangle(a_.second.position, b_.position, c_.position).normal();
						if (dot(tn, mn) > 0)
							pm.add_triangle(pmp::Vertex(a), pmp::Vertex(b), pmp::Vertex(c));
						else
							pm.add_triangle(pmp::Vertex(a), pmp::Vertex(c), pmp::Vertex(b));
					}
				}
			}
		}

		// run the regularization
		pmp::SurfaceRemeshing rms(pm);
		rms.uniform_remeshing(1, 10);

		{
			NavigationData n;
			// vertices
			for (const auto &v : pm.positions())
			{
				NavigationPoint p;
				p.position = vec3(v[0], v[1], v[2]);
				n[numeric_cast<uint32>(n.size())] = templates::move(p);
			}
			// edges
			for (const auto &e : pm.edges())
			{
				uint32 a = pm.vertex(e, 0).idx();
				uint32 b = pm.vertex(e, 1).idx();
				n[a].neighbors.insert(b);
				n[b].neighbors.insert(a);
			}
			std::swap(n, navigation);
		}
	}
}

void optimizeNavigation(NavigationData &navigationParam, const NavmeshOptimizationConfig &cfg)
{
	std::swap(navigation, navigationParam);
	validateNavigationDebug();
	joinIdenticalVertices();
	printNavigationStatistics();
	prepareOrigData();
	pmpRegularization();
	updateVertexPropertiesFromOrig();
	markBorderVertices();
	printNavigationStatistics();
	for (uint32 iteration = 0; iteration < cfg.iterations; iteration++)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "navmesh optimizing iteration: " + iteration);
		optimizeVertexPositions();
		addNeighborEdges();
		splitLongEdges();
		joinCloseVertices();
		removeSuboptimalEdges();
		removeEmptyVertices();
		updateVertexPropertiesFromOrig();
		printNavigationStatistics();
	}
	clearOrigData();
	validateNavigationRelease();
	std::swap(navigation, navigationParam);
}
