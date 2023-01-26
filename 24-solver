from itertools import product, permutations

# todo: generalize to more digits

ops = [None]*4
ops[0] = ["{0}{1}{2}{3}"]
ops[1] = ["{0}{4}{1}{2}{3}",
          "{0}{1}{4}{2}{3}",
          "{0}{1}{2}{4}{3}"]
ops[2] = ["({0}{1}{4}{2}){5}{3}",
          "{0}{1}{4}({2}{5}{3})",
          "({0}{4}{1}{2}){5}{3}",
          "{0}{4}({1}{2}{5}{3})",
          "({0}{4}{1}){5}{2}{3}",
          "{0}{4}({1}{5}{2}{3})"]
ops[3] = ["(({0}{4}{1}){5}{2}){6}{3}",
          "({0}{4}({1}{5}{2})){6}{3}",
          "({0}{4}{1}){5}({2}{6}{3})",
          "{0}{4}({1}{5}({2}{6}{3}))",
          "{0}{4}(({1}{5}{2}){6}{3})"]

def find(digs,oper,res):
    for i in range(4):
        for op in ops[i]:
            for pdigs in permutations(digs):
                for poper in product(oper,repeat=i):
                    x = op.format(*pdigs,*poper)
                    try:
                        if abs(eval(x) - res) < 1e-9:
                            print(x)
                    except:
                        pass

find([1,4,5,9],["+","-","*","/"],4)
