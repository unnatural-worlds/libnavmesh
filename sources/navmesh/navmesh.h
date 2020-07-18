#include <unnatural-navmesh/navmesh.h>

#include <set>
#include <map>
#include <vector>

using namespace unnatural;

struct NavigationPoint
{
	std::set<uint32> neighbors;
	vec3 position;
	vec3 normal;
	uint8 terrain = 0; // terrain type index
	bool border = false;
};

typedef std::map<uint32, NavigationPoint> NavigationData;

void optimizeNavigation(NavigationData &navigation, const NavmeshOptimizationConfig &cfg);
