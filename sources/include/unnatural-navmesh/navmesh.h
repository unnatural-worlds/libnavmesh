#include <cage-core/polyhedron.h>

namespace unnatural
{
	using namespace cage;

	struct NavmeshOptimizationConfig
	{
		float targetScale = 1;
		uint32 iterations = 20;
		bool pmpRegularization = true;
		bool markBorderVertices = true;
	};

	Holder<Polyhedron> navmeshOptimize(const Holder<Polyhedron> &nav, const NavmeshOptimizationConfig &cfg);
}
