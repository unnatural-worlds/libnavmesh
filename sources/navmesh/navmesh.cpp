#include "navmesh.h"

#include <cage-core/meshAlgorithms.h>

namespace unnatural
{
	using namespace navoptim;

	Holder<Mesh> navmeshOptimize(const Holder<Mesh> &navigation, const NavmeshOptimizeConfig &config)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "initiating navigation optimization");

		switch (navigation->type())
		{
		case MeshTypeEnum::Triangles:
		case MeshTypeEnum::Lines:
			break;
		default:
			CAGE_THROW_ERROR(Exception, "unsupported mesh type");
		}

		Holder<Mesh> nav = navigation->copy();
		
		{
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "merging close vertices");
			MeshMergeCloseVerticesConfig cfg;
			cfg.distanceThreshold = 1e-3;
			meshMergeCloseVertices(+nav, cfg);
		}

		Graph graph = convertMeshToGraph(+nav, config.tileSize);
		graphValidationUnconditional(graph); // input validation
		printStatistics(graph);
		const SpatialGraph original = convertToSpatialGraph(graph);

		if (nav->type() == MeshTypeEnum::Triangles && config.pmpRegularization)
		{
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "initial regularization");
			MeshRegularizeConfig cfg;
			cfg.targetEdgeLength = config.tileSize;
			cfg.iterations = 10;
			meshRegularize(+nav, cfg);
			graph = convertMeshToGraph(+nav, config.tileSize);
			graphValidationDebugOnly(graph);
			printStatistics(graph);
		}
		else
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "skipping initial regularization");

		if (nav->type() == MeshTypeEnum::Triangles && config.markBorderVertices)
			markBorderNodes(graph);
		else
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "skipping marking border nodes");

		for (uint32 iteration = 0; iteration < config.iterations; iteration++)
		{
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", Stringizer() + "navmesh optimizing iteration: " + iteration);
			optimizeNodePositions(graph);
			addNeighborEdges(graph);
			splitLongEdges(graph);
			joinCloseNodes(graph);
			removeSuboptimalEdges(graph);
			updateNodeProperties(graph, original);
			printStatistics(graph);
		}

		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "finishing navigation optimization");
		graphValidationUnconditional(graph); // output validation
		nav = convertGraphToMesh(graph, config.tileSize);
		return nav;
	}
}

