#ifndef guard_navmesh_private_h_sd5rf4gzhdrtg
#define guard_navmesh_private_h_sd5rf4gzhdrtg

#include <vector>

#include <cage-core/enumerate.h>
#include <cage-core/flatSet.h>
#include <cage-core/geometry.h>
#include <cage-core/mesh.h>
#include <cage-core/spatialStructure.h>
#include <unnatural-navmesh/navmesh.h>

namespace navoptim
{
	using namespace unnatural;

	struct Node
	{
		Vec3 position;
		Vec3 normal;
		uint8 terrain = 0; // terrain type index
		bool border = false;
	};

	struct Graph
	{
		std::vector<Node> nodes;
		std::vector<FlatSet<uint32>> neighbors;
	};

	struct SpatialGraph : public Graph
	{
		Holder<SpatialStructure> spatialData;
		Holder<SpatialQuery> spatialQuery;
	};

	Graph convertMeshToGraph(const Mesh *poly, Real tileSize);
	Holder<Mesh> convertGraphToMesh(const Graph &graph, Real tileSize);
	const SpatialGraph convertToSpatialGraph(const Graph &graph);

	void printStatistics(const Graph &graph);
	void graphValidationUnconditional(const Graph &graph);
	void graphValidationDebugOnly(const Graph &graph);

	std::vector<uint32> sharedNeighbors(const Graph &graph, uint32 a, uint32 b);
	FlatSet<uint32> connectedNodes(const Graph &graph, const FlatSet<uint32> &start);
	uint32 totalEdgesCount(const Graph &graph);
	void removeEmptyNodes(Graph &graph);

	void markBorderNodes(Graph &graph);

	void optimizeNodePositions(Graph &graph, Real factor);
	void snapNodesToCollider(Graph &graph, const Collider *collider, Real tileSize);
	void addNeighborEdges(Graph &graph);
	void splitLongEdges(Graph &graph);
	void joinCloseNodes(Graph &graph);
	void removeSuboptimalEdges(Graph &graph);
	void updateNodeProperties(Graph &graph, const SpatialGraph &original);
}

#endif // guard_navmesh_private_h_sd5rf4gzhdrtg
