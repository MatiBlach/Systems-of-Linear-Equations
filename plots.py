import matplotlib.pyplot as plt
import pandas as pd

#Zadanie B
data = pd.read_csv('dataB.csv', header=None)
jacobi_data = data.iloc[0].values
gauss_data = data.iloc[1].values
plt.plot(jacobi_data, label='Jacobi')
plt.plot(gauss_data, label='Gauss-Seidel')
plt.xlabel('Iterations')
plt.ylabel('Residual Norm')
plt.title('Change in Residual Norm over Iterations')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.savefig('wykresB.png')
plt.show()

#Zadanie C
data = pd.read_csv('dataC.csv', header=None)
jacobi_data = data.iloc[0].values
gauss_data = data.iloc[1].values
plt.plot(jacobi_data, label='Jacobi')
plt.plot(gauss_data, label='Gauss-Seidel')
plt.xlabel('Iterations')
plt.ylabel('Residual Norm')
plt.title('Change in Residual Norm over Iterations')
plt.yscale('log')  
plt.grid(True)
plt.legend()
plt.savefig('wykresC.png')
plt.show()

#ZadanieE
data = pd.read_csv('dataE.csv')

# Separate data for each method
n_sizes = data['n_size']
jm_time = data['Jacobi_Method_Time']
gsm_time = data['Gauss_Seidel_Method_Time']
lu_time = data['LU_Factorization_Time']

# Create plot
plt.plot(n_sizes, jm_time, label='Jacobi Method')
plt.plot(n_sizes, gsm_time, label='Gauss-Seidel Method')
plt.plot(n_sizes, lu_time, label='LU Factorization')
plt.xlabel('Matrix Size')
plt.ylabel('Time (seconds)')
plt.title('Time vs. Size')
plt.legend()
plt.grid(True)
plt.savefig('wykresE.png')
plt.show()