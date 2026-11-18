from lib import SEQ
print(len(SEQ))
for i in range(0, len(SEQ), 10):
    print(i + 1, SEQ[i:i + 10])
