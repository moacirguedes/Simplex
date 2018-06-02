using System;

namespace Simplex
{
    class Gauss
    {
        public double[] CalcGauss(int size, double[,] mat)
        {
            double[] result = new double[10];
            double c, som = 0;
            for (int j = 0; j < size; j++)
            {
                for (int i = 0; i < size; i++)
                {
                    if (i > j)
                    {
                        c = mat[i, j] / mat[j, j];
                        for (int k = 0; k <= size; k++)
                        {
                            mat[i, k] = mat[i, k] - c * mat[j, k];
                        }
                    }
                }
            }
            result[size - 1] = mat[size - 1, size] / mat[size - 1, size - 1];

            for (int i = size - 2; i >= 0; i--)
            {
                som = 0;
                for (int j = i; j < size; j++)
                {
                    som = som + mat[i, j] * result[j];
                }
                result[i] = (mat[i, size] - som) / mat[i, i];
            }
            return result;
        }
    }
}
