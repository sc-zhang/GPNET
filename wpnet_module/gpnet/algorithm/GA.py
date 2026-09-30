from numpy import argmax, array, concatenate, inf, random, zeros, arange
from gpnet.io.data_io import DataSaver


class GA:
    def __init__(
        self,
        type_info,
        allele_name,
        outfile,
        func,
        is_lower_better,
        is_store_iter,
        seed=None,
        **kwargs,
    ):
        self.__type_info = array(type_info, dtype=int)
        self.__allele_name = allele_name
        self.__outfile = outfile + ".iter"
        self.__func = func
        self.__is_lower_better = is_lower_better
        self.__is_store_iter = is_store_iter
        self.__n_generations = kwargs.get("n_generations", 200)
        self.__pop_size = kwargs.get("pop_size", 100)
        self.__cross_rate = kwargs.get("cross_rate", 0.8)
        self.__init_mutation_rate = kwargs.get("mutation_rate", 0.05)
        self.__mutation_rate = self.__init_mutation_rate
        self.__n_features = len(type_info)

        if seed:
            random.seed(seed)
        else:
            random.seed()

        self.data = None
        self.__pop = None
        self.__iter_data = []
        self.__iter_score = []

    def __to_onehot(self, individual):
        onehot_segments = []

        for idx, cat_count in zip(individual, self.__type_info):
            segment = zeros(cat_count, dtype=int)
            segment[idx] = 1
            onehot_segments.append(segment)

        return concatenate(onehot_segments)

    def __generate_new(self, pop, fitness):
        # tournament selection
        tournament = random.randint(0, self.__pop_size, size=(self.__pop_size, 3))
        winner = tournament[
            arange(self.__pop_size), argmax(fitness[tournament], axis=1)
        ]
        selected = pop[winner]

        # elitism
        elite_idx = argmax(fitness)
        elite = pop[elite_idx].copy()

        new_pop = [elite]

        # crossover and mutation
        while len(new_pop) < self.__pop_size:
            a, b = random.randint(0, self.__pop_size, size=2)
            parent_a = selected[a]
            parent_b = selected[b]
            child = parent_a.copy()
            # crossover
            if random.rand() < self.__cross_rate:
                mask = random.random(self.__n_features) < 0.5
                child[mask] = parent_b[mask]
            # mutation
            for i, n_cat in enumerate(self.__type_info):
                if random.rand() < self.__mutation_rate:
                    old = child[i]
                    if n_cat > 1:
                        choices = arange(n_cat)
                        choices = choices[choices != old]
                        child[i] = random.choice(choices)
            new_pop.append(child)
        return array(new_pop)

    def run(self):

        self.__pop = array(
            [random.randint(0, self.__type_info) for _ in range(self.__pop_size)]
        )

        best_fitness = -inf
        best_individual = None

        for _ in range(self.__n_generations):

            scores = array([self.__func(self.__to_onehot(ind)) for ind in self.__pop])
            fitness = -scores if self.__is_lower_better else scores
            current = argmax(fitness)

            if fitness[current] > best_fitness:
                best_fitness = fitness[current]
                best_individual = self.__pop[current].copy()
                best_score = scores[current]
                self.__iter_data.append(self.__to_onehot(best_individual))
                self.__iter_score.append(best_score)
            self.__pop = self.__generate_new(self.__pop, fitness)
            self.__mutation_rate = max(
                0.02, self.__init_mutation_rate * (1 - _ * 1.0 / self.__n_generations)
            )

        self.data = self.__to_onehot(best_individual)

        if self.__is_store_iter:
            ds = DataSaver(self.__outfile)
            ds.save_iter(
                self.__type_info,
                self.__allele_name,
                self.__iter_data,
                self.__iter_score,
            )
