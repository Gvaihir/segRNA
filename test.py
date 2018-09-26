import sys

# print usage message
__doc__ %= sys.argv[0]

with open(sys.argv[1], 'r') as file:
    contents = file.read()
print(contents)