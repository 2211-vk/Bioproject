import flwr as fl  # pip install flwr
from deap import base, creator, tools
import numpy as np
import random

# GA setup (local client)
def create():
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

def evaluate(individual, data):  # Fitness: Predict resistance from path
    score = sum(individual) * np.mean(data['gc_content'])  # Simplified
    return score,

def toolbox_setup(df_meta):
    create()
    toolbox = base.Toolbox()
    toolbox.register("individual", tools.initRepeat, creator.Individual, random.random, n=10)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", evaluate, data=df_meta)
    return toolbox

# Flower Client
class GAFlowerClient(fl.client.NumPyClient):
    def __init__(self, toolbox, df_meta):
        self.toolbox = toolbox_setup(df_meta)

    def get_parameters(self, config):
        return [random.random() for _ in range(10)]  # GA params as weights

    def fit(self, parameters, config):
        # Local GA evolution
        population = [creator.Individual(random.random() for _ in range(10)) for _ in range(50)]
        for gen in range(5):  # 5 generations
            offspring = tools.selTournament(population, len(population), 2)
            offspring = [self.toolbox.clone(ind) for ind in offspring]
            for ind in offspring:
                self.toolbox.evaluate(ind)
            population[:] = offspring
        return [ind.fitness.values[0] for ind in population[:10]], 100, {}  # Return top params

# Start server
fl.server.start_server(strategy=fl.server.strategy.FedAvg(), config=fl.server.ServerConfig(num_rounds=3))
