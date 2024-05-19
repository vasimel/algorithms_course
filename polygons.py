import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from scipy.spatial import distance
import networkx as nx

# Функция для проверки, пересекается ли ячейка с полигоном
def is_cell_intersecting_polygon(cell_center, cell_size, polygon):
    cell_polygon = Polygon([
        (cell_center[0] - cell_size / 2, cell_center[1] - cell_size / 2),
        (cell_center[0] + cell_size / 2, cell_center[1] - cell_size / 2),
        (cell_center[0] + cell_size / 2, cell_center[1] + cell_size / 2),
        (cell_center[0] - cell_size / 2, cell_center[1] + cell_size / 2)
    ])
    return cell_polygon.intersects(polygon)

# Задаем полигоны и точки начала и конца
polygons = [Polygon([(1, 0), (3, 0), (2.5, 2.5), (0.6, 2)]),
        Polygon([(3, 3), (4, 3), (4, 4), (3, 4)])]
start_point = Point(0, 0)
end_point = Point(5, 5)

# Определяем границы пространства
min_x = min(min(polygon.bounds[0] for polygon in polygons), start_point.x, end_point.x)
max_x = max(max(polygon.bounds[2] for polygon in polygons), start_point.x, end_point.x)
min_y = min(min(polygon.bounds[1] for polygon in polygons), start_point.y, end_point.y)
max_y = max(max(polygon.bounds[3] for polygon in polygons), start_point.y, end_point.y)

# Создаем сетку ячеек
cell_size = 0.5
x_coords = np.arange(min_x, max_x, cell_size)
y_coords = np.arange(min_y, max_y, cell_size)

# Создаем граф
graph = nx.Graph()

# Добавляем узлы и ребра в граф
for x in x_coords:
    for y in y_coords:
        cell_center = (x + cell_size / 2, y + cell_size / 2)
        if not any(is_cell_intersecting_polygon(cell_center, cell_size, polygon) for polygon in polygons):
            graph.add_node(cell_center)
            if x > min_x:
                left_neighbor = (x - cell_size / 2, y + cell_size / 2)
                if left_neighbor in graph:
                    graph.add_edge(cell_center, left_neighbor)
            if y > min_y:
                bottom_neighbor = (x + cell_size / 2, y - cell_size / 2)
                if bottom_neighbor in graph:
                    graph.add_edge(cell_center, bottom_neighbor)

# Добавляем начальную и конечную точки в граф
graph.add_node(start_point)
graph.add_node(end_point)
for node in graph.nodes:
    if node != start_point and node != end_point:
        if distance.euclidean(node, (start_point.x, start_point.y)) < cell_size * np.sqrt(2):
            graph.add_edge(start_point, node)
        if distance.euclidean(node, (end_point.x, end_point.y)) < cell_size * np.sqrt(2):
            graph.add_edge(end_point, node)


# Ищем кратчайший путь
path = nx.shortest_path(graph, start_point, end_point)

# Визуализация
plt.figure(figsize=(8, 8))
for polygon in polygons:
    plt.plot(*polygon.exterior.xy, color='red')
# Изменяем способ передачи координат точек в plt.plot
path_coords = [(p.x, p.y) if isinstance(p, Point) else p for p in path]
plt.plot(*zip(*path_coords), color='blue', marker='o')

plt.xlim(min_x, max_x)
plt.ylim(min_y, max_y)
plt.grid(True)
plt.show()

