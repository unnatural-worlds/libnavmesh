
#include <cage-core/polyhedron.h>

namespace unnatural
{
	using namespace cage;

	struct NavmeshOptimizationConfig
	{
		PolyhedronRegularizationConfig regularization;

		NavmeshOptimizationConfig();
	};

	Holder<Polyhedron> navmeshOptimize(Holder<Polyhedron> &&nav, const NavmeshOptimizationConfig &cfg);
}
