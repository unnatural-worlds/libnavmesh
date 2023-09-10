#include <set>

#include "navmesh.h"

namespace navoptim
{
	Graph convertMeshToGraph(const Mesh *poly, Real tileSize)
	{
		Graph res;

		{ // vertices
			const Real scale = 1 / tileSize;
			const bool hasUv = !poly->uvs().empty();
			const uint32 vc = poly->verticesCount();
			res.nodes.reserve(vc);
			for (uint32 i = 0; i < vc; i++)
			{
				Node n;
				n.position = poly->position(i) * scale;
				n.normal = poly->normal(i);
				if (hasUv)
				{
					const Vec2 uv = saturate(poly->uv(i));
					n.terrain = numeric_cast<uint8>(uv[0] * 32);
					n.border = uv[1] > 0.5;
				}
				res.nodes.push_back(n);
			}
			res.neighbors.resize(vc);
			for (auto &it : res.neighbors)
				it.reserve(10);
		}

		// indices
		switch (poly->type())
		{
			case MeshTypeEnum::Triangles:
			{
				const uint32 tc = poly->indicesCount() / 3;
				const auto inds = poly->indices();
				for (uint32 i = 0; i < tc; i++)
				{
					const uint32 a = inds[i * 3 + 0];
					const uint32 b = inds[i * 3 + 1];
					const uint32 c = inds[i * 3 + 2];
					auto &na = res.neighbors[a];
					auto &nb = res.neighbors[b];
					auto &nc = res.neighbors[c];
					na.insert(b);
					na.insert(c);
					nb.insert(a);
					nb.insert(c);
					nc.insert(a);
					nc.insert(b);
				}
			}
			break;
			case MeshTypeEnum::Lines:
			{
				const uint32 lc = poly->indicesCount() / 2;
				const auto inds = poly->indices();
				for (uint32 i = 0; i < lc; i++)
				{
					const uint32 a = inds[i * 2 + 0];
					const uint32 b = inds[i * 2 + 1];
					res.neighbors[a].insert(b);
					res.neighbors[b].insert(a);
				}
			}
			break;
			default:
				CAGE_THROW_ERROR(Exception, "unsupported mesh type");
		}

		return res;
	}

	Holder<Mesh> convertGraphToMesh(const Graph &graph, Real tileSize)
	{
		Holder<Mesh> res = newMesh();
		res->type(MeshTypeEnum::Lines);

		{ // positions
			std::vector<Vec3> v;
			v.reserve(graph.nodes.size());
			for (const auto &it : graph.nodes)
				v.push_back(it.position * tileSize);
			res->positions(v);
		}

		{ // normals
			std::vector<Vec3> v;
			v.reserve(graph.nodes.size());
			for (const auto &it : graph.nodes)
				v.push_back(it.normal);
			res->normals(v);
		}

		{ // uvs
			std::vector<Vec2> v;
			v.reserve(graph.nodes.size());
			for (const auto &it : graph.nodes)
			{
				Vec2 u;
				u[0] = (Real(it.terrain) + 0.5) / 32;
				u[1] = it.border;
				v.push_back(u);
			}
			res->uvs(v);
		}

		{ // indices
			std::vector<uint32> inds;
			inds.reserve(graph.nodes.size() * 10);
			for (const auto &it : enumerate(graph.neighbors))
			{
				const uint32 a = numeric_cast<uint32>(it.index);
				for (const uint32 b : *it)
				{
					if (b > a)
					{
						inds.push_back(a);
						inds.push_back(b);
					}
				}
			}
			res->indices(inds);
		}

		return res;
	}

	const SpatialGraph convertToSpatialGraph(const Graph &graph)
	{
		SpatialGraph res;
		res.nodes = graph.nodes;
		res.neighbors = graph.neighbors;
		res.spatialData = newSpatialStructure({});
		res.spatialQuery = newSpatialQuery(res.spatialData.share());
		for (const auto &it : enumerate(res.nodes))
			res.spatialData->update(numeric_cast<uint32>(it.index), Aabb(it->position));
		res.spatialData->rebuild();
		return res;
	}

	void printStatistics(const Graph &graph)
	{
		Real lenSum;
		Real lenMin = Real::Infinity();
		Real lenMax;
		uint32 nghCnt = 0;
		for (const auto &a : enumerate(graph.neighbors))
		{
			for (uint32 b : *a)
			{
				const Real d = distance(graph.nodes[a.index].position, graph.nodes[b].position);
				lenSum += d;
				lenMin = min(lenMin, d);
				lenMax = max(lenMax, d);
				nghCnt++;
			}
		}
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "edge length minimum: " + lenMin + ", maximum: " + lenMax);
		Real lenAvg = lenSum / nghCnt;
		Real devSum1 = 0;
		Real devSum2 = 0;
		for (const auto &a : enumerate(graph.neighbors))
		{
			for (uint32 b : *a)
			{
				const Real d = distance(graph.nodes[a.index].position, graph.nodes[b].position);
				devSum1 += abs(d - lenAvg);
				devSum2 += abs(d - 1);
			}
		}
		Real devAvg1 = devSum1 / nghCnt;
		Real devAvg2 = devSum2 / nghCnt;
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "edge length average: " + lenAvg + ", deviation: " + devAvg1 + ", unfitness: " + devAvg2);
	}

	void graphValidationUnconditional(const Graph &graph)
	{
		if (graph.nodes.empty())
			CAGE_THROW_ERROR(Exception, "validation error: graph is empty");
		if (graph.neighbors.size() != graph.nodes.size())
			CAGE_THROW_ERROR(Exception, "validation error: inconsistent array sizes");
		for (const auto &n : graph.nodes)
		{
			if (!n.position.valid())
				CAGE_THROW_ERROR(Exception, "validation error: invalid node position");
			if (!n.normal.valid() || abs(lengthSquared(n.normal) - 1) > 1e-3)
				CAGE_THROW_ERROR(Exception, "validation error: invalid node normal");
		}
		for (const auto &a : enumerate(graph.neighbors))
		{
			if (a->empty())
				CAGE_THROW_ERROR(Exception, "validation error: node has no edges");
			for (uint32 b : *a)
			{
				if (b >= graph.nodes.size())
					CAGE_THROW_ERROR(Exception, "validation error: invalid edge index");
				if (graph.neighbors[b].count(numeric_cast<uint32>(a.index)) != 1)
					CAGE_THROW_ERROR(Exception, "validation error: edge is missing counterpart");
				if (a.index == b)
					CAGE_THROW_ERROR(Exception, "validation error: edge to itself");
			}
		}
		{ // check that the whole graph is single connected component
			std::vector<bool> visited;
			visited.resize(graph.nodes.size(), false);
			visited[0] = true;
			std::set<uint32> open;
			open.insert(0);
			while (!open.empty())
			{
				const uint32 n = *open.begin();
				open.erase(open.begin());
				CAGE_ASSERT(visited[n]);
				for (uint32 i : graph.neighbors[n])
				{
					if (visited[i])
						continue;
					visited[i] = true;
					open.insert(i);
				}
			}
			for (bool v : visited)
			{
				if (!v)
					CAGE_THROW_ERROR(Exception, "validation error: disconnected component");
			}
		}
	}

	void graphValidationDebugOnly(const Graph &graph)
	{
#ifdef CAGE_DEBUG
		graphValidationUnconditional(graph);
#endif // CAGE_DEBUG
	}
}
