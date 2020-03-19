"""
Uses unittest and compares this intersection code to Carel Faber's MgnLb.txt maple program.

To use, you'll need to  edit this file so that fabers_dir is the correct value.

Then, there are two choices:

1)  Edit the gn_pairs_to_check attribute to include the pairs you want to check.  Then, run
$ sage testvsMgnLb.py
and it will check them all.

2) Create a file in the format:

g, n, [list of classes by Faber's index]

e.g.:

1,3, [1,2,3]
0,4, [1]

and then do
$ sage testvsMgnLb.py testlist
(where testlist is the filename) and this will just check the things in your list.  This will make a quick check that you can run often.
"""
fabers_dir = "/home/drew/Dropbox/fabers maple code/Send10/"
gn_pairs_tocheck = ( (1,1), (0,4), (0,5), (0,6), (1,2), (1,3), (1,4), (2,0), (2,1), (2,2), (3,0))#( (5,0), (6,0), (7,0) )
    #(4,0)
    #( (1,1), (0,4), (0,5), (0,6), (1,2), (1,3), (1,4), (2,0), (2,1), (2,2), (3,0))



import  unittest
from sage.all import maple
from intersection_ui import intnum
from TautRing3 import Mgn
from get_all_combis import get_all_valid_tuples


class TestVsMaple(unittest.TestCase):

    def check(self, g, n, c, verbose = True):
        #M = Mgn(g,n)
        #print M
        #C = Combinations( [M.fabers_classes_indexes()] * M.dimension, M.dimension )
        #print C
        #for c in Combinations( M.fabers_classes_indexes() * M.dimension, M.dimension ):
            #c = c[0] #????
        if verbose:
            print("Checking on M_{0}, {1}, classes {2}...".format(g, n, c))
        myans = intnum(g,n,c, confirm = False)
        #print "mgn({0}, {1})".format(g, c)
        mapleans = maple("mgn({0}, {1}, {2})".format(g, n, c))
        if verbose:
            print("values {3}, {4}".format(g,n, c, myans, mapleans))
        self.assertEqual(myans, mapleans, "FAILED:  for genus {0}, {1} marked points with classes {2}, my anwers: {3},  maple answer: {4}".format(g,n, c,myans,mapleans) )





    def start_Maple(self):
        maple('currentdir("{0}")'.format(fabers_dir))
        maple.eval('read "MgnLb.txt"')
        print("Loaded MgnLb.")

    def test_from_file(self, filename):
        self.start_Maple()
        with open(filename) as f:
            for line in f:
                g,n,c = eval(line)
                self.check(g,n,c)

    def test_all(self):
        self.start_Maple()
        for g,n in gn_pairs_tocheck:
            print("Checking {0}, {1}...".format(g, n))
            M = Mgn(g,n)
            for c in get_all_valid_tuples(M.degree_index_dict(), M.dimension):
                self.check(g,n,c)
        #save_data("genus4567test")

if __name__ == '__main__':
    #global filename
    import sys
    #tester = TestVsMaple()
    if len(sys.argv) > 1:
        #TestVsMaple.filename = sys.argv[1]
        tester = TestVsMaple("test_from_file")
        tester.test_from_file(sys.argv[1])
    else:
        tester = TestVsMaple("test_all")
        tester.test_all()
