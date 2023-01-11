x = 10
y = x

if id(x) == id(y):
	print("x and y refer to the same object")




x = 10
y = x
x += 1

if id(x) != id(y):
	print("x and y do not refer to the same object")



def func():
		
	# All these variables get memory
	# allocated on stack
	a = 20
	b = []
	c = ""

i = 0

while i < 100:
    i = i + 1clear
    
