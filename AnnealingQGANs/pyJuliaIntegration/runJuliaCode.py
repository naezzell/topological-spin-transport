from julia.api import Julia
import time
from julia import Main


jpath = "/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia"
jl = Julia(runtime=jpath, compiled_modules=False)

jl.eval('include("juliaForLoop.jl")')

n = int(1e8)
Main.n = n

tic = time.time()
jlstr = f"for_loop(n)"
result = jl.eval(jlstr)
toc = time.time()
juliaTime = toc - tic
print(f"Julia for loop time: {juliaTime}\n")

def for_loop(iterations = 100):
   a = 0
   for i in range(iterations):
       a = a+1
   return a

tic = time.time()
result = for_loop(n)
toc = time.time()
pythonTime = toc - tic
print(f"Python for loop time: {pythonTime}\n")
