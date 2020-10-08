#include <cage-core/polyhedron.h>

namespace unnatural
{
	using namespace cage;

	struct NavmeshOptimizationConfig
	{
		float tileSize = 10; // the goal, by default, is that average distance between any two neighbor tiles is 10 meters
		uint32 iterations = 20;
		bool pmpRegularization = true;
		bool markBorderVertices = true;
	};

	Holder<Polyhedron> navmeshOptimize(const Holder<Polyhedron> &nav, const NavmeshOptimizationConfig &cfg);
}
