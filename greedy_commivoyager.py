# решение задачи коммивояжера жадным алгоритмом. 
# Граф представлен матрицей смежности



def find_nearest_neighbor(graph, current_city, visited):
    min_distance = float('inf')
    nearest_neighbor = None
    for neighbor, distance in enumerate(graph[current_city]):
        if neighbor not in visited and distance < min_distance:
            min_distance = distance
            nearest_neighbor = neighbor
    return nearest_neighbor

def greedy_tsp(graph, start_city=0):
    num_cities = len(graph)
    visited = {start_city}
    path = [start_city]
    total_distance = 0

    current_city = start_city
    while len(visited) < num_cities:
        next_city = find_nearest_neighbor(matrix, current_city, visited)
        visited.add(next_city)
        path.append(next_city)
        total_distance += matrix[current_city][next_city]
        current_city = next_city

    # из последнего узла возвращаемся в начальный
    total_distance += matrix[current_city][start_city]
    path.append(start_city)

    return path, total_distance

# пример
matrix = [
    [0, 40, 15, 20],
    [40, 0, 35, 25],
    [15, 35, 0, 30],
    [20, 25, 30, 0]
]

path, total_distance = greedy_tsp(matrix)
print(f"Path: {path}")
print(f"Total Distance: {total_distance}")
