# pylint: disable=all
import blueshift

test_data = "../../tests/DCinput1.csv"
test_info = "../../tests/DCinfo1.csv"

############
# Example 1 - Supervised Internal Standard + Nearest Neighbor Drift Correction
############

b = blueshift.DriftCorrection(test_data, test_info)

# output full raw dataset
print(b.data)

# output reference pool file names
print(b.pools)

# perform internal standard correction on rdata using the given internal standard(s)
b.internal_standard_correct(internal_standards="15R-15-methyl PGA2")

# output full drift corrected dataset
print(b.data)

# calculate CVs for IS corrected data
b.calc_cvs()

# output CVs for IS corrected data
print(b.cvs)

# perform nearest neighbor correction on cdata using the given pool
b.pool_correct(interpolation="NN", pool="PREFB", override=False)

# output full drift corrected dataset
print(b.data)

# calculate CVs for all metabolites
b.calculate_cvs()

# save the corrected data and CVs as a csv with the provided filename
b.to_csv("IS-NN_correct.csv")
b.report()

############
# Example 2 - Supervised Linear Drift Correction
############

b = blueshift.DriftCorrection(test_data, test_info)

# output full raw dataset
print(b.data)

# performs the indicated smooth method on cdata (equivalent to rdata until one drift correction method is done)
b.pool_correct(interpolation="linear", override=False, pool="PREFA")

# calculate CVs for all metabolites
b.calculate_cvs()

# save the corrected data and CVs as a csv with a default filename
b.to_csv()
b.report()

############
# Example 3 - Unsupervised Internal Standard + Linear Drift Correction with Override
############

b = blueshift.DriftCorrection(test_data, test_info)

# automatically perform IS correction followed by linear correction with overrides
b.correct(
    internal_standard_kwargs={"internal_standards": "15R-15-methyl PGA2"},
    pool_kwargs={"interpolation": "smooth", "override": True, "pool": "PREFA"},
)

# calculate CVs for all metabolites
b.calculate_cvs()

# output full drift corrected dataset
print(b.data)

# save the corrected data and CVs as a csv with a default filename
b.to_csv()
b.report()
