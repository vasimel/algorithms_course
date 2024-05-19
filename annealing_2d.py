import numpy as np

# Задаем функцию и ее производную
def f(x):
    return x**4 + x**3 - 1

def grad_f(x):
    return 4*x**3 + 3*x**2

# Функция для градиентного спуска
def gradient_descent(x, learning_rate, max_iter, x_min=-10, x_max=10):
    for _ in range(max_iter):
        x = x - learning_rate * grad_f(x)
        x = np.clip(x, x_min, x_max)  # Ограничиваем x
    return x

# Функция для метода отжига
def simulated_annealing(initial_temp, cooling_rate, max_iter):
    current_temp = initial_temp
    x = np.random.uniform(-10, 10)  # Начальная точка в диапазоне от -10 до 10
    current_solution = x
    best_solution = x
    best_energy = f(x)

    for i in range(max_iter):
        new_solution = current_solution + np.random.normal(0, 1)
        new_energy = f(new_solution)


        # сохраняем лучшее решение!!
        if new_energy < best_energy:
            best_solution = new_solution
            best_energy = new_energy

        # если новое решение лучше предыдущего, прыгаем в него,
        # если нет, то мы в него прыгаем в какой-то вероятностью, которая с уменьшением температуры уменьшается
        if new_energy < f(current_solution) or np.exp((f(current_solution) - new_energy) / current_temp) > np.random.rand():
            current_solution = new_solution

        current_temp *= cooling_rate

        # Применяем градиентный спуск на каждом шаге
        current_solution = gradient_descent(current_solution, 0.01, 10)

    return best_solution, best_energy

# Параметры метода отжига
initial_temp = 1000
cooling_rate = 0.99
max_iter = 1000

# Запускаем метод отжига
best_solution, best_energy = simulated_annealing(initial_temp, cooling_rate, max_iter)

print(f"Лучшее найденное решение: {best_solution}")
print(f"Значение функции в лучшем решении: {best_energy}")
