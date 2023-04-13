print("How many rows? : ")
n = int(input())

xs = []
ys = []

print("Copy the table and paste : ")

for i in range(n):
    aa = list(map(float,input().split())) 
    a = aa[0]
    b = aa[1]
    xs.append(a)
    ys.append(b)

print(xs)
print(ys)