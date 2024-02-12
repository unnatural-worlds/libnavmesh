#ifndef guard_navmesh_h_awsegersd4gh
#define guard_navmesh_h_awsegersd4gh

#include <cage-core/collider.h>
#include <cage-core/mesh.h>

namespace unnatural
{
	using namespace cage;

	struct NavmeshOptimizeConfig
	{
		const Mesh *navigation = nullptr;
		const Collider *collider = nullptr;
		Real tileSize = 10; // the goal, by default, is that average distance between any two neighbor tiles is 10 meters
		uint32 iterations = 20;
		bool pmpRegularization = true;
		bool markBorderVertices = true;
	};

	Holder<Mesh> navmeshOptimize(NavmeshOptimizeConfig &config);
}

#endif // guard_navmesh_h_awsegersd4gh
