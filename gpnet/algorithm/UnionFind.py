class UnionFind():
    def __init__(self, size):
        self.__f = [i for i in range(0, size)]

    def find(self, x):
        if self.__f[x] == x:
            return x
        self.__f[x] = self.find(self.__f[x])
        return self.__f[x]

    def union(self, x, y):
        fx = self.find(x)
        fy = self.find(y)
        if fx != fy:
            self.__f[fy] = fx
