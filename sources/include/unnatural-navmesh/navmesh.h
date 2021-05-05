#ifndef guard_navmesh_h_awsegersd4gh
#define guard_navmesh_h_awsegersd4gh

#include <cage-core/core.h>

namespace unnatural
{
	using namespace cage;

	struct NavmeshOptimizeConfig
	{
		float tileSize = 10; // the goal, by default, is that average distance between any two neighbor tiles is 10 meters
		uint32 iterations = 20;
		bool pmpRegularization = true;
		bool markBorderVertices = true;
	};

	Holder<Mesh> navmeshOptimize(const Holder<Mesh> &navigation, const NavmeshOptimizeConfig &config);
}

#endif // guard_navmesh_h_awsegersd4gh
