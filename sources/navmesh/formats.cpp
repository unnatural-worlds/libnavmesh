#include "navmesh.h"

namespace navoptim
{
	Graph convertMeshToGraph(const Polyhedron *poly, real tileSize)
	{
		Graph res;

		{ // vertices
			const real scale = 1 / tileSize;
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
					n.terrain = numeric_cast<uint8>(poly->uv(i)[0] * 32);
					n.border = poly->uv(i)[1] > 0.5;
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
		case PolyhedronTypeEnum::Triangles:
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
		} break;
		case PolyhedronTypeEnum::Lines:
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
		} break;
		default:
			CAGE_THROW_ERROR(Exception, "unsupported polyhedron type");
		}

		return res;
	}

	Holder<Polyhedron> convertGraphToMesh(const Graph &graph, real tileSize)
	{
		Holder<Polyhedron> res = newPolyhedron();
		res->type(PolyhedronTypeEnum::Lines);

		{ // positions
			std::vector<vec3> v;
			v.reserve(graph.nodes.size());
			for (const auto &it : graph.nodes)
				v.push_back(it.position * tileSize);
			res->positions(v);
		}

		{ // normals
			std::vector<vec3> v;
			v.reserve(graph.nodes.size());
			for (const auto &it : graph.nodes)
				v.push_back(it.normal);
			res->normals(v);
		}

		{ // uvs
			std::vector<vec2> v;
			v.reserve(graph.nodes.size());
			for (const auto &it : graph.nodes)
			{
				vec2 u;
				u[0] = (real(it.terrain) + 0.5) / 32;
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
		res.spatialQuery = newSpatialQuery(+res.spatialData);
		for (const auto &it : enumerate(res.nodes))
			res.spatialData->update(numeric_cast<uint32>(it.index), aabb(it->position));
		res.spatialData->rebuild();
		return res;
	}

	void printStatistics(const Graph &graph)
	{
		real lenSum;
		real lenMin = real::Infinity();
		real lenMax;
		uint32 nghCnt = 0;
		for (const auto &a : enumerate(graph.neighbors))
		{
			for (uint32 b : *a)
			{
				const real d = distance(graph.nodes[a.index].position, graph.nodes[b].position);
				lenSum += d;
				lenMin = min(lenMin, d);
				lenMax = max(lenMax, d);
				nghCnt++;
			}
		}
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "edge length minimum: " + lenMin + ", maximum: " + lenMax);
		real lenAvg = lenSum / nghCnt;
		real devSum1 = 0;
		real devSum2 = 0;
		for (const auto &a : enumerate(graph.neighbors))
		{
			for (uint32 b : *a)
			{
				const real d = distance(graph.nodes[a.index].position, graph.nodes[b].position);
				devSum1 += abs(d - lenAvg);
				devSum2 += abs(d - 1);
			}
		}
		real devAvg1 = devSum1 / nghCnt;
		real devAvg2 = devSum2 / nghCnt;
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "edge length average: " + lenAvg + ", deviation: " + devAvg1 + ", unfitness: " + devAvg2);
	}

	void graphValidationUnconditional(const Graph &graph)
	{
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
	}

	void graphValidationDebugOnly(const Graph &graph)
	{
#ifdef CAGE_DEBUG
		graphValidationUnconditional(graph);
#endif // CAGE_DEBUG
	}
}
