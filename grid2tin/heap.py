# Binary heap implementation based on Sedgewick and Wayne's 
# IndexMaxPQ.java: http://algs4.cs.princeton.edu/24pq/IndexMaxPQ.java

import array


class Heap:
    def __init__(self):
        self.pq = array.array('i', [-1])
        self.qp = array.array('i', [-1])
        self.keys = array.array('f', [-1])
        self.elements = [None]
        self.N = 0
        self.I = 0

    def insert(self, key, element):
        self.N += 1
        self.I += 1
        self.pq.append(self.I)
        self.qp.append(self.N)
        self.keys.append(key)
        self.elements.append(element)

        self.swim(self.N)

        return self.I

    def max(self):
        return self.keys[self.pq[1]], self.elements[self.pq[1]]

    def del_max(self):
        maximum = self.pq[1]
        self.exchange(1, self.N)
        self.pq.pop()
        self.N -= 1

        self.sink(1)
        self.qp[maximum] = -1

    def pop(self):
        maximum = self.max()
        self.del_max()
        return maximum

    def less(self, i, j):
        return self.keys[self.pq[i]] < self.keys[self.pq[j]]

    def swim(self, k):
        while k > 1 and self.less(k // 2, k):
            self.exchange(k, k // 2)
            k //= 2

    def sink(self, k):
        while 2 * k <= self.N:
            j = 2 * k
            if j < self.N and self.less(j, j + 1):
                j += 1
            if not self.less(k, j):
                break
            self.exchange(k, j)
            k = j

    def exchange(self, i, j):
        temp = self.pq[i]
        self.pq[i] = self.pq[j]
        self.pq[j] = temp

        self.qp[self.pq[i]] = i
        self.qp[self.pq[j]] = j

    def delete(self, i):
        index = self.qp[i]
        self.exchange(index, self.N)
        self.pq.pop()
        self.N -= 1
        if index <= self.N:
            self.swim(index)
            self.sink(index)
        self.elements[i] = None
        self.qp[i] = -1
