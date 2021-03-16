#include "navmesh.h"

namespace unnatural
{
	using namespace navoptim;

	Holder<Polyhedron> navmeshOptimize(const Holder<Polyhedron> &navigation, const NavmeshOptimizationConfig &config)
	{
		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "initiating navigation optimization");

		switch (navigation->type())
		{
		case PolyhedronTypeEnum::Triangles:
		case PolyhedronTypeEnum::Lines:
			break;
		default:
			CAGE_THROW_ERROR(Exception, "unsupported polyhedron type");
		}

		Holder<Polyhedron> nav = navigation->copy();
		
		{
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "merging close vertices");
			PolyhedronCloseVerticesMergingConfig cfg;
			cfg.distanceThreshold = 1e-3;
			polyhedronMergeCloseVertices(+nav, cfg);
		}

		Graph graph = convertMeshToGraph(+nav, config.tileSize);
		printStatistics(graph);
		const SpatialGraph original = convertToSpatialGraph(graph);

		if (nav->type() == PolyhedronTypeEnum::Triangles && config.pmpRegularization)
		{
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "initial regularization");
			PolyhedronRegularizationConfig cfg;
			cfg.targetEdgeLength = config.tileSize;
			cfg.iterations = 10;
			polyhedronRegularize(+nav, cfg);
			graph = convertMeshToGraph(+nav, config.tileSize);
			printStatistics(graph);
		}
		else
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "skipping initial regularization");

		if (nav->type() == PolyhedronTypeEnum::Triangles && config.markBorderVertices)
			markBorderNodes(graph);
		else
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", "skipping marking border nodes");

		for (uint32 iteration = 0; iteration < config.iterations; iteration++)
		{
			CAGE_LOG(SeverityEnum::Info, "libnavmesh", stringizer() + "navmesh optimizing iteration: " + iteration);
			optimizeNodePositions(graph);
			addNeighborEdges(graph);
			splitLongEdges(graph);
			joinCloseNodes(graph);
			removeSuboptimalEdges(graph);
			updateNodeProperties(graph, original);
			printStatistics(graph);
		}

		CAGE_LOG(SeverityEnum::Info, "libnavmesh", "finishing navigation optimization");
		graphValidationUnconditional(graph);
		nav = convertGraphToMesh(graph, config.tileSize);
		return nav;
	}
}

