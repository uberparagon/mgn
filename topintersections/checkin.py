"""
This is very useful for debugging recursive programs.  It is not a part of the intersection theory code.
"""


level = 0
def checkin(func):
    def checkin_func(*args):
        global level
        print("." * level + "Entering {0} with arguments {1}".format(func.__name__, args))
        level += 1
        value = func(*args)
        level -= 1
        print("." * level + "Exiting {0} with arguments {1}, returning value {2}".format(func.__name__, args, value))

        return value
    return checkin_func






class Node:
    def __init__(self, data, parent = None, func = "unknown"):
        self.children = []
        self.data =  str(data)
        self.real_data = data
        self.parent = parent
        self.value = None
        self.func = func
        self.flag = ""

    def add_child(self, data, func="unknown"):
        child = Node(data, self,func)
        self.children.append(child)
        return child

    def make_big_maple_code(self, filename):
        with open(filename, "w") as f:
            f.write("read MgnF;")
            for c in self.children:
                f.write(c.maple_code() + "\n")

    def maple_code(self):
        M = self.data[0][0]
        m = self.data[1]
        expon_list = []
        for l in ([M.num(d)]*expon for d, expon in m.decompose_monomial()):
            expon_list += l
        return "mgn(" + str(M.genus) + ", "+ str( expon_list )  + ");\n"

    def __str__(self):
        _str = repr(self)
        for c,i in zip(self.children, range(len(self.children))):
            _str += "\n[{0}] ".format(i) + repr(c)
        return _str

    def inFaberNotation(self): #not done yet!
        _str = repr(self)
        for c,i in zip(self.children, range(len(self.children))):
            _str += "\n[{0}] ".format(i) + repr(c)
        return _str

    def __repr__(self):
        return self.func + " " + str(self.data) + " = " + str(self.value) + " " + self.flag

    def __getitem__(self,i):
        return self.children[i]

    def mark_ancestors(self, message):
        if self.parent != None:
            self.parent.flag += ", " +  message
            self.parent.mark_ancestors(message)


def cireset():
    global root, current, level
    root = Node("root")
    current = root
    level = 0
cireset()

def breadth_first_tree(func):
    def bft(*args):
        global root, current, level
        #print "." * level + "Entering {0} with arguments {1}".format(func.__name__, args)
        level += 1
        current = current.add_child(args, func.__name__)
        value = func(*args)
        current.value = value
        current = current.parent
        level -= 1
        #print "." * level + "Exiting {0} with arguments {1}, returning value {2}".format(func.__name__, args, value)
        return value

    return bft

fabers_dir = "/home/drew/Dropbox/fabers maple code/Send10/"

from sage.all import maple
maple('currentdir("{0}")'.format(fabers_dir))
maple.eval('read "MgnLb.txt"')
print("Loaded MgnLb.")

def breadth_first_tree_MgnLb(func):
    def bft(*args):
        global root, current, level
        #print "." * level + "Entering {0} with arguments {1}".format(func.__name__, args)
        level += 1
        current = current.add_child(args, func.__name__)

        g= args[0].space.genus
        n= args[0].space.n
        c= args[0].get_MgnLb_indexes()
        mapleans = maple("mgn({0}, {1}, {2})".format(g, n, c))
        value = func(*args)
        if mapleans != value:
            current.flag = "wrong... should be " + str(c) + " = " + str(mapleans)
            current.mark_ancestors("bad child")
        current.value = value
        current = current.parent
        level -= 1
        #print "." * level + "Exiting {0} with arguments {1}, returning value {2}".format(func.__name__, args, value)
        return value

    return bft
