using System;

namespace Com_Methods
{
    /// <summary>
    /// преобразования Гивенса
    /// </summary>
    class Givens_Transformation
    {
        /// <summary>
        /// реализация процедуры ортогонализации по методу вращений Гивенса
        /// </summary>
        /// <param name="A - исходная матрица"></param>
        /// <param name="Q - ортогональная матрица преобразований"></param>
        /// <param name="R - результат"></param>
        public static void Givens_Orthogonalization(Matrix A, Matrix Q, Matrix R)
        {
            double help1, help2;

            //косинус, синус
            double c = 0, s = 0;

            //запись матрицы А в R
            R.Copy(A);

            //алгоритм вращения Гивенса: для каждого столбца
            for (int j = 0; j < R.N - 1; j++)
            {
                //просматриваем строки в столбце
                for (int i = j + 1; i < R.M; i++)
                {
                    //если очередной элемент под диагональю не нулевой, то требуется поворот вектора
                    if (Math.Abs(R.Elem[i][j]) > CONST.EPS)
                    {
                        help1 = Math.Sqrt(Math.Pow(R.Elem[i][j], 2) + Math.Pow(R.Elem[j][j], 2));
                        c = R.Elem[j][j] / help1;
                        s = R.Elem[i][j] / help1;

                        //A_new = Gt * A
                        for (int k = j; k < R.N; k++)
                        {
                            help1 = c * R.Elem[j][k] + s * R.Elem[i][k];
                            help2 = c * R.Elem[i][k] - s * R.Elem[j][k];
                            R.Elem[j][k] = help1;
                            R.Elem[i][k] = help2;
                        }

                        //перемножаем строки матрицы Q на трансп.матрицу преобразования Q = Q * G
                        for (int k = 0; k < Q.M; k++)
                        {
                            help1 = c * Q.Elem[k][j] + s * Q.Elem[k][i];
                            help2 = c * Q.Elem[k][i] - s * Q.Elem[k][j];
                            Q.Elem[k][j] = help1;
                            Q.Elem[k][i] = help2;
                        }
                    }
                }
            }
        }
    }
}
