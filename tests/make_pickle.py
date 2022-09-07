from bmxp.eclipse import MSAligner
import pickle

a = MSAligner("test1.csv", "test2.csv", names=["HP1", "HP2"])
a.align()
file = open("mseclipse.pickle", "wb")
pickle.dump(a, file)
file.close()

