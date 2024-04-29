#calculating the age of a sample using the Carbon-14 isotope knowing that the half life of carbon-14 is 5700 years
    

 # Half-life of Carbon-14 in years
import math
half_life_c14 = 5700

    # Assuming remaining ratio of Carbon-14 in the sample
remaining_ratio = float(input("Enter the remaining ratio of Carbon-14 in the sample (between 0 and 1): "))

# Calculate the age using the formula: age = -(half_life) * ln(remaining_ratio) / ln(2)
age = -(half_life_c14) * math.log(remaining_ratio) / math.log(2)
print (age)

print(f"The age of the sample is approximately {age:.2f} years.")
    


