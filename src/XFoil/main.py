import subprocess

XFOIL_PATH = "bin/XFoil/XFoil.exe"
Re = 1000000
Alpha = 0

# Create a new XFoil process
xfoil = subprocess.Popen(XFOIL_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

xfoil.stdin.write("NACA 2412\n")
xfoil.stdin.write("OPER\n")
xfoil.stdin.write(f"Visc {Re}\n")
xfoil.stdin.write(f"ALFA {Alpha}\n")
xfoil.stdin.write("DUMP data/BL.dat\n")
xfoil.stdin.write("QUIT\n")

stdout, stderr = xfoil.communicate()

print("=====================")
print(stdout)
print("=====================")
print(stderr)
print("=====================")
