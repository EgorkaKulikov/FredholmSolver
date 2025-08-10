using System;
using System.Collections.Generic;

using MathNet.Numerics.Integration;

namespace FredholmSolver
{
  //Алиасы типов для упрощения кода
  using Matrix = Dictionary<Tuple<int, int>, double>;
  using Vector = Dictionary<int, double>;

  public static class MathUtils
  {
    /// <summary>
    /// Обёртка для решателя системы линейных уравнений.
    /// 
    /// Так как у нас все индексы лежат в интервале -2..N-1,
    /// а C# не поддерживает отрицательные индексы, приходится 
    /// сдвигать их все на 2, решать систему, а потом делать обратный сдвиг.
    /// </summary>
    public static Vector SolveSystem(Matrix matrix, Vector vector)
    {
      double[,] arrayMatrix = new double[Configuration.GridPoints + 2, Configuration.GridPoints + 2];
      double[] arrayVector = new double[Configuration.GridPoints + 2];
      var coeffs = new Dictionary<int, double>();

      //сдвиг индексов вперёд
      foreach (var index in matrix.Keys)
      {
        arrayMatrix[index.Item1 + 2, index.Item2 + 2] = matrix[index];
      }
      foreach (int index in vector.Keys)
      {
        arrayVector[index + 2] = vector[index];
      }

      //решение системы
      var solution = SolveSystem(arrayMatrix, arrayVector);

      //сдвиг индексов назад
      for (int i = 0; i < solution.Length; i++)
      {
        coeffs[i - 2] = solution[i];
      }

      return coeffs;
    }

    /// <summary>
    /// Численное интегрирование функции func на отрезке [left, right]
    /// Используется библиотека MathNet.Numeric
    /// </summary>
    public static double Integrate(Func<double, double> func, double left, double right)
    {
      var pointsCount = 1000;
      var ncIntegral = NewtonCotesTrapeziumRule.IntegrateComposite(func, left, right, pointsCount);
      var smIntegral = SimpsonRule.IntegrateComposite(func, left, right, pointsCount);
      return 0.5 * (ncIntegral + smIntegral);
    }

    /// <summary>
    /// Решение системы линейных уравнений методом Гаусса
    /// </summary>
    private static double[] SolveSystem(double[,] matrix, double[] vector)
    {
      double[] x = new double[vector.Length];

      for (int k = 0; k < vector.Length - 1; k++)
      {
        for (int i = k + 1; i < vector.Length; i++)
        {
          for (int j = k + 1; j < vector.Length; j++)
          {
            matrix[i, j] = matrix[i, j] - matrix[k, j] * (matrix[i, k] / matrix[k, k]);
          }
          vector[i] = vector[i] - vector[k] * matrix[i, k] / matrix[k, k];
        }
      }

      for (int k = vector.Length - 1; k >= 0; k--)
      {
        double sum = 0.0;
        for (int j = k + 1; j < vector.Length; j++)
        {
          sum = sum + matrix[k, j] * x[j];
        }

        x[k] = (vector[k] - sum) / matrix[k, k];
      }

      return x;
    }
  }
}
