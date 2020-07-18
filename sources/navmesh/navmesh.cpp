#include "unnatural-navmesh/navmesh.h"

namespace unnatural
{
	NavmeshOptimizationConfig::NavmeshOptimizationConfig()
	{

	}

	Holder<Polyhedron> navmeshOptimize(Holder<Polyhedron> &&nav, const NavmeshOptimizationConfig &cfg)
	{
		nav->regularize(cfg.regularization);
		return templates::move(nav);
	}
}

