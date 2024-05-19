// даны непересекающиеся полигоны. 
//Снова плоскость. Заданы непересекающиеся полигоны, точка начала и целевая точка. 
//Требуется построить оптимальный (в смысле построения) маршрут из точки в точку. 
//Расстояние между точками - просто L2 расстояние. В первую очередь каким-то способом строится 
//граф, а затем при помощи алгоритма Дийкстры или поиска в глубину уже находится оптимальный 
//маршрут.
//Если мы разобьем все пространство на ячейки, то центры таких ячеек - будут узлы графа. Если ячейка пересекается с каким-то полигоном, мы считаем, что в ячейку нельзя заходить
//- ячеечный алгоритм для преобразования непрерывного ограниченного полигонами пространства в граф – 4 балла



#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

struct Point {
    double x, y;

    Point() : x(0), y(0) {} // Конструктор по умолчанию
    Point(double x, double y) : x(x), y(y) {}

    double distance(const Point& other) const {
        return std::sqrt(std::pow(x - other.x, 2) + std::pow(y - other.y, 2));
    }
};

struct Cell {
    Point center;
    bool isBlocked;

    Cell(double x, double y) : center(x, y), isBlocked(false) {}
};

bool isCellIntersectingPolygon(const Cell& cell, const std::vector<Point>& polygon) {

    for (size_t i = 0; i < polygon.size(); i++) {
        if (polygon[i].x == cell.center.x && polygon[i].y == cell.center.y) {
            return true;
        }
    }

    // Проверка на пересечение ячейки с полигоном. Если центр ячейки лежит внутри полигона, 
    // то ячейка пересекается с полигоном. Тут метод пересекающего луча.
    int count = 0;
    for (size_t i = 0; i < polygon.size(); i++) {
        
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % polygon.size()];
        // проверка, лежат ли концы ребра полигона по разные стороны от точки. 
        if ((p1.y <= cell.center.y && cell.center.y < p2.y) ||
            (p2.y <= cell.center.y && cell.center.y < p1.y)) {
                // находим точку пересечения луча из точки в направлении оси x с этим ребром полигона
                // по подобию треугольников
            double x = p1.x + (cell.center.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
            if (x == cell.center.x) {
                return true;
            }
            else if (x > cell.center.x) {
                count++;
            }
        }
    }
    return count % 2 != 0; // Нечётное количество пересечений означает, что точка внутри полигона
}

using namespace boost;
    typedef adjacency_list<vecS, vecS, undirectedS, Point, property<edge_weight_t, double>> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef std::pair<int, int> Edge;

Vertex find_vertex_by_point(const Point& point, const Graph& graph, const std::vector<Vertex>& vertexMap) {
    for (const auto& v : vertexMap) {
        if (graph[v].x == point.x && graph[v].y == point.y) {
            return v;
        }
    }
    return boost::graph_traits<Graph>::null_vertex();
}


int main() {


    // Определение полигонов и точек начала и конца
    std::vector<std::vector<Point>> polygons = {
        {{1, 0}, {3, 0}, {2.5, 2.5}, {0.6, 2}},
        {{3, 3}, {4, 3}, {4, 4}, {3, 4}}
    };
    Point startPoint(0, 0);
    Point endPoint(5, 5);

    // Определение границ пространства
    double minX = std::min({startPoint.x, endPoint.x});
    double maxX = std::max({startPoint.x, endPoint.x});
    double minY = std::min({startPoint.y, endPoint.y});
    double maxY = std::max({startPoint.y, endPoint.y});
    for (const auto& polygon : polygons) {
        for (const auto& point : polygon) {
            minX = std::min(minX, point.x);
            maxX = std::max(maxX, point.x);
            minY = std::min(minY, point.y);
            maxY = std::max(maxY, point.y);
        }
    }

    // Создание сетки
    double cellSize = 0.5;
    std::vector<std::vector<Cell>> grid;
    for (double x = minX - cellSize; x < maxX; x += cellSize) {
        std::vector<Cell> row;
        for (double y = minY - cellSize; y < maxY; y += cellSize) {
            row.emplace_back(x + cellSize / 2, y + cellSize / 2);
        }
        grid.push_back(row);
    }

    // Отметка заблокированных ячеек
    for (auto& row : grid) {
        for (auto& cell : row) {
            for (const auto& polygon : polygons) {
                if (isCellIntersectingPolygon(cell, polygon)) {
                    cell.isBlocked = true;
                    break;
                }
            }
        }
    }


    for (auto& row : grid) {
        for (auto& cell : row) {
            std::cout << cell.center.x << ", " << cell.center.y << ", " << cell.isBlocked << std::endl;
        }
    }
    // Создание графа, коллекции с вершинами, ребрами и весами
    Graph graph;
    std::vector<Vertex> vertexMap;
    std::vector<Edge> edges;
    std::vector<double> weights;


size_t vertexIndex = 0;
for (size_t i = 0; i < grid.size(); ++i) {
    for (size_t j = 0; j < grid[i].size(); ++j) {
        // заблоченные ячейки добавляем в граф, но их не связываем ни с чем
        Vertex v = add_vertex(grid[i][j].center, graph);
        vertexMap.push_back(v);
            // связываем соседние ячейки в одном столбце (нижнюю с верхней)
        if (i > 0 && !grid[i - 1][j].isBlocked && !grid[i][j].isBlocked) {
            edges.emplace_back(vertexMap[vertexIndex], vertexMap[vertexIndex - grid[i].size()]);
            weights.push_back(cellSize);
        }
            // связываем соседние ячейки в одном ряду (правую с левой)
        if (j > 0 && !grid[i][j - 1].isBlocked && !grid[i][j].isBlocked) {
            edges.emplace_back(vertexMap[vertexIndex], vertexMap[vertexIndex - 1]);
            weights.push_back(cellSize);
        }
        vertexIndex++;
    }
}


 // создаем соответствующие ребра уже в графе
    for (size_t i = 0; i < edges.size(); ++i) {
        add_edge(edges[i].first, edges[i].second, weights[i], graph);
    }



// Попытка найти начальную и конечную точки в графе
Vertex startVertex = find_vertex_by_point(startPoint, graph, vertexMap);
if (startVertex == boost::graph_traits<Graph>::null_vertex()) {
    std::cout << "Start vertex not found: creating one" << std::endl;
    startVertex = add_vertex(startPoint, graph);
}

Vertex endVertex = find_vertex_by_point(endPoint, graph, vertexMap);
if (endVertex == boost::graph_traits<Graph>::null_vertex()) {
    std::cout << "End vertex not found: creating one" << std::endl;
    endVertex = add_vertex(endPoint, graph);
}

// Соединение начальной и конечной точек с ближайшими вершинами графа
for (const auto& v : vertexMap) {
    Point p = graph[v];
    if (startPoint.distance(p) < cellSize * std::sqrt(2) && v != startVertex) {
        add_edge(startVertex, v, startPoint.distance(p), graph);
    }
    if (endPoint.distance(p) < cellSize * std::sqrt(2) && v != endVertex) {
        add_edge(endVertex, v, endPoint.distance(p), graph);
    }
}


for (const auto &v : vertexMap) { 
    std::cout << graph[v].x << " " << graph[v].y << std::endl;
}
    std::vector<Vertex> predecessors(num_vertices(graph));
    std::vector<double> distances(num_vertices(graph));
    dijkstra_shortest_paths(graph, startVertex,
                            predecessor_map(&predecessors[0]).distance_map(&distances[0]));

if (predecessors[endVertex] != endVertex) { // Проверяем, существует ли путь
    // Извлечение пути
    std::vector<Point> path;
    for (Vertex v = endVertex; v != startVertex; v = predecessors[v]) {
        path.push_back(graph[v]);
    }
    path.push_back(startPoint);

    // Print the path
    std::cout << "Shortest path:" << std::endl;
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        std::cout << "(" << it->x << ", " << it->y << ")" << std::endl;
    }
} else {
    std::cout << "No path exists between the start and end points." << std::endl;
}

    return 0;
}
