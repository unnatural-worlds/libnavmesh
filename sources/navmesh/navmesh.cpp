#include "navmesh.h"

namespace unnatural
{
	Holder<Polyhedron> navmeshOptimize(const Holder<Polyhedron> &nav, const NavmeshOptimizationConfig &cfg_)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "navigation optimization");

		NavmeshOptimizationConfig cfg = cfg_;
		NavigationData data;

		{ // load
			const uint32 vc = nav->verticesCount();

			{ // vertices
				const bool hasUv = !nav->uvs().empty();
				const real scale = 1 / cfg.tileSize;
				for (uint32 i = 0; i < vc; i++)
				{
					auto &v = data[i];
					v.position = nav->position(i) * scale;
					v.normal = nav->normal(i);
					if (hasUv)
					{
						v.terrain = numeric_cast<uint8>(nav->uv(i)[0] * 32);
						v.border = nav->uv(i)[1] > 0.5;
					}
				}
			}

			// indices
			switch (nav->type())
			{
			case PolyhedronTypeEnum::Triangles:
			{
				const uint32 tc = nav->indicesCount() / 3;
				const auto inds = nav->indices();
				for (uint32 i = 0; i < tc; i++)
				{
					const uint32 a = inds[i * 3 + 0];
					const uint32 b = inds[i * 3 + 1];
					const uint32 c = inds[i * 3 + 2];
					auto &na = data[a].neighbors;
					auto &nb = data[b].neighbors;
					auto &nc = data[c].neighbors;
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
				const uint32 lc = nav->indicesCount() / 2;
				const auto inds = nav->indices();
				for (uint32 i = 0; i < lc; i++)
				{
					const uint32 a = inds[i * 2 + 0];
					const uint32 b = inds[i * 2 + 1];
					auto &na = data[a].neighbors;
					auto &nb = data[b].neighbors;
					na.insert(b);
					nb.insert(a);
				}
				cfg.pmpRegularization = false;
				cfg.markBorderVertices = false;
			} break;
			default:
				CAGE_THROW_ERROR(Exception, "unsupported polyhedron type");
			}
		}

		optimizeNavigation(data, cfg);

		{ // save
			Holder<Polyhedron> res = newPolyhedron();
			res->type(PolyhedronTypeEnum::Lines);

			{ // positions
				std::vector<vec3> v;
				v.reserve(data.size());
				for (const auto &it : data)
					v.push_back(it.second.position * cfg.tileSize);
				res->positions(v);
			}

			{ // normals
				std::vector<vec3> v;
				v.reserve(data.size());
				for (const auto &it : data)
					v.push_back(it.second.normal);
				res->normals(v);
			}

			{ // uvs
				std::vector<vec2> v;
				v.reserve(data.size());
				for (const auto &it : data)
				{
					vec2 u;
					u[0] = (real(it.second.terrain) + 0.5) / 32;
					u[1] = it.second.border;
					v.push_back(u);
				}
				res->uvs(v);
			}

			{ // indices
				std::vector<uint32> inds;
				inds.reserve(data.size());
				for (const auto &it : data)
				{
					const uint32 a = it.first;
					for (const uint32 b : it.second.neighbors)
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
	}
}

