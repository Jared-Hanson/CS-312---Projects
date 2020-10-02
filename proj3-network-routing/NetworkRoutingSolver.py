#!/usr/bin/python3


from CS312Graph import *
import time


# Simple class to store information along with reference to a specific node.
# Will be used in the hash-dictionaries and to store distances in the queues
class arrayObject(object):
    def __init__(self, node):
        self.dist = float("inf")
        self.prev = None
        self.node = node
        self.index = -1



#  Unsorted list implementation.
class Unsorted_List(object):

    # Uses two data structures:
    # An unsorted dictionary (hashmap) of key pairs [Node_id] -> ArrayObject containing specific node
    # An unsorted dictionary (nodes) of key paris [NodeID] -> Distance. This is the main priority queue
    def __init__(self, nodeList):
        self.hashmap = {}
        self.nodes = {}
        for index, node in enumerate(nodeList):
            self.hashmap[node.node_id] = arrayObject(node)

    # Fill queue with references and keys
    def makeQueue(self, sourceNodeId):
        for x in self.hashmap:
            self.nodes[x] = float('inf')
        self.nodes[sourceNodeId] = 0

    # Finds lowest distance and deletes/returns it.
    def deletemin(self):
        minId = -1
        minDist = -1
        for x in self.nodes:
            if minId == -1:
                minId = x
                minDist = self.nodes[x]
            elif self.nodes[x] < minDist:
                minId = x
                minDist = self.nodes[x]
        # Update dist value in hashmap
        tempNode = self.hashmap[minId]
        tempNode.dist = minDist
        self.nodes.pop(minId)
        return tempNode

    # Inserts new node into end of array
    def insert(self, node):
        self.hashmap[node.node_id] = arrayObject(node)
        self.nodes[node.node_id] = float('inf')

    # Acceses node through the hashmap and lowers its distance value.
    def decreaseKey(self, curNode_ID, distance, nodeV_ID):
        self.nodes[nodeV_ID] = distance
        self.hashmap[nodeV_ID].dist = distance
        self.hashmap[nodeV_ID].prev = self.hashmap[curNode_ID]


# Second Implementation of Priority Queue, Uses a Bianary Treee/Heap
class bianaryHeap(object):

    # Uses two data structures:
    # An unsorted dictionary (hashmap) of key pairs [Node_id] -> ArrayObject containing specific node
    # An unsorted list (nodes) of nodes. This is the main priority queue
    def __init__(self, nodeList):
        self.nodes = []
        self.hashmap = {}
        for index, node in enumerate(nodeList):
            self.hashmap[node.node_id] = arrayObject(node)

    # initializes both lists, puts the source node at the front of the main queue. Fills all others in a line behind it.
    def makeQueue(self, sourceNodeId):
        self.nodes.append(self.hashmap[sourceNodeId])
        self.nodes[0].dist = 0
        self.nodes[0].index = 0
        for x in self.hashmap:
            if x != sourceNodeId:
                self.nodes.append(self.hashmap[x])
                self.nodes[len(self.nodes) - 1].index = len(self.nodes) - 1

    # Given a node and an index, the function rises that node up the tree swapping nodes until the tree is corectly balanced
    def bubbleup(self, nodeX, i):
        i = i + 1
        p = i // 2
        while i > 1 and self.nodes[p - 1].dist > nodeX.dist:
            self.nodes[i - 1] = self.nodes[p - 1]
            self.nodes[i - 1].index = i - 1
            i = p
            p = i//2
        self.nodes[i - 1] = nodeX
        self.nodes[i - 1].index = i - 1

    # Uses bubble up function to insert a node a the end of the list and bubble it up into its correct place
    def insert(self, nodeX):
        self.nodes.append(nodeX)
        # I store the index of the specific node to increase speed when locating the node with a hashmap
        # This keeps the hashmap up to date with the queue.
        self.nodes[len(self.nodes) - 1].index = len(self.nodes) - 1
        self.bubbleup(nodeX, len(self.nodes) - 1)

    # Locates the given node through the hashmap, lowers its value, and bubbles it into place
    def decreaseKey(self, curNode_ID, distance, nodeV_ID):
        nodeX = self.hashmap[nodeV_ID]
        nodeX.dist = distance
        nodeX.prev = self.hashmap[curNode_ID]
        self.bubbleup(nodeX, nodeX.index)

    # Takes the lowest node (index 0) and returns it. It then swaps the last node into the first
    # place and sifts it down.
    def deletemin(self):
        if len(self.nodes) == 0:
            return None
        else:
            x = self.nodes[0]
            self.siftdown(self.nodes.pop(), 0)
            return x

    # Looks at nodes smallest child with (minchild) and swaps it with the node if it is less.
    # continues until the tree is balanced
    def siftdown(self, nodeX, i):
        i = i + 1
        c = self.minChild(i)
        while c != 0 and self.nodes[c-1].dist < nodeX.dist:
            self.nodes[i - 1] = self.nodes[c - 1]
            self.nodes[i - 1].index = i - 1
            i = c
            c = self.minChild(i)

        # Takes care of boundry case where the heap only has one entry
        if len(self.nodes) > 1:
            self.nodes[i - 1] = nodeX
            self.nodes[i - 1].index = i - 1

    # Looks at the two children of the node and returns a 0 if it has no children less than it
    # else, returns the index of the min child
    def minChild(self, i):
        if 2 * i > len(self.nodes):
            return 0
        elif i * 2 + 1 > len(self.nodes):
            return i * 2
        else:
            if self.nodes[(i * 2) - 1].dist < self.nodes[(i * 2 + 1) - 1].dist:
                return i * 2
            else:
                return i * 2 + 1








class NetworkRoutingSolver:
    def __init__( self ):
        self.queue = None
        pass

    def initializeNetwork( self, network ):
        assert( type(network) == CS312Graph )
        self.network = network
    # Uses the linked list of previous nodes to create lines on the graph
    def getShortestPath( self, destIndex ):
        self.dest = destIndex

        path_edges = []
        curNode = self.queue.hashmap[destIndex]
        prevNode = curNode.prev
        distance = float('inf')
        # Continues tracing previous nodes and creating lines until source node is reached
        while prevNode is not None:
            path_edges.append( (curNode.node.loc, prevNode.node.loc, '{:.0f}'.format(curNode.dist - prevNode.dist)) )
            curNode = prevNode
            prevNode = curNode.prev
            distance = self.queue.hashmap[destIndex].dist
        return {'cost':distance, 'path':path_edges}
    # Initialize heap based on boolean indicator and run Dikstras algoryythm.
    def computeShortestPaths( self, srcIndex, use_heap=False ):
        self.source = srcIndex
        t1 = time.time()
        # Select correct bianaryHeap
        if use_heap:
            self.queue = bianaryHeap(self.network.nodes)
            self.queue.makeQueue(srcIndex)
        else:
            self.queue = Unsorted_List(self.network.nodes)
            self.queue.makeQueue(srcIndex)

        # Dikstras Algorythm
        iteration = 0
        while len(self.queue.nodes) != 0:
            u = self.queue.deletemin()
            for edge in u.node.neighbors:
                if edge.length + u.dist < self.queue.hashmap[edge.dest.node_id].dist:
                    self.queue.hashmap[edge.dest.node_id].dist = u.dist + edge.length
                    self.queue.decreaseKey(u.node.node_id, self.queue.hashmap[edge.dest.node_id].dist, edge.dest.node_id)
            iteration += 1

        t2 = time.time()
        return (t2-t1)


