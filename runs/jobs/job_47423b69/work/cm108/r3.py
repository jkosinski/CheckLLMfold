from lib import *
import r2

REMOVE = {(5,20),(4,19),(3,20),(2,15),(33,60),(34,59),(35,59),(76,99),(56,70),(37,95),(34,95),(13,58),(13,66),
          (49,106),(26,41),(26,70),(29,33),(1,4),(10,13),(29,38),(28,34)}
ADD = {(5,45),(5,49),(7,42),(5,46),            # helix 34-49 hydrophobic side reaching I5/L7
       (27,99),(25,102),(23,106),(23,103),     # helix 95-107 hydrophobic side reaching I23/V25/F27
       (30,75),(22,54),(82,107),(9,65),(12,69),(45,54),   # close in the round-2 model and sensible
       (17,54),(16,54),(17,24)}                # V16/L17 on the T54/L24 side of the strands

def contacts():
    s = set(norm(r2.contacts())) - set(norm(REMOVE))
    return s | set(norm(ADD))

if __name__ == '__main__':
    run(contacts(), 'r3')
