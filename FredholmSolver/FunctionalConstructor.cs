using System;
using System.Collections.Generic;

namespace FredholmSolver
{
  public class FunctionalConstructor
  {
    /// <summary>
    /// Вычисление значения аппроксимационного функционала.
    /// 
    /// Значение определяется типом и номером выбранного функционала, 
    /// функцией, к которой он примененён, а также основной и дополнительной сеткой.
    /// </summary>
    public double FunctionalValue(
      int funcNumber,
      Func<double, double> f,
      Dictionary<int, double> grid,
      Dictionary<int, double> additionalGrid
      )
    {
      int j = funcNumber;
      switch (Configuration.ApproxType)
      {
        case ApproximationType.ShoenbergMarsden:
          return f(additionalGrid[j]);

        case ApproximationType.Sablonniere:
          if (j == -2 || j == Configuration.GridPoints - 1)
          {
            return f(additionalGrid[j]);
          }

          double a = -1.0 / 8;
          double b = 5.0 / 4;
          double c = -1.0 / 8;
          return a * f(additionalGrid[j - 1]) + b * f(additionalGrid[j]) + c * f(additionalGrid[j + 1]);

        case ApproximationType.Trigonometric:
          if (j == -2 || j == Configuration.GridPoints - 1)
          {
            return f(additionalGrid[j]);
          }

          double S1 = 1.0 * (Math.Cos(grid[j + 1]) - Math.Cos(grid[j + 2])) / Math.Sin(grid[j + 2] - grid[j + 1]);
          double S2 = 1.0 * (Math.Sin(grid[j + 2]) - Math.Sin(grid[j + 1])) / Math.Sin(grid[j + 2] - grid[j + 1]);

          double zn = (Math.Sin(additionalGrid[j + 1]) - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j]) - Math.Cos(additionalGrid[j - 1]))
            - (Math.Sin(additionalGrid[j]) - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j + 1]) - Math.Cos(additionalGrid[j - 1]));

          double ch_b = (Math.Sin(additionalGrid[j + 1]) - Math.Sin(additionalGrid[j - 1])) * (S2 - Math.Cos(additionalGrid[j - 1]))
            - (S1 - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j + 1]) - Math.Cos(additionalGrid[j - 1]));

          double ch_c = (S1 - Math.Sin(additionalGrid[j - 1])) * (Math.Cos(additionalGrid[j]) - Math.Cos(additionalGrid[j - 1]))
            - (Math.Sin(additionalGrid[j]) - Math.Sin(additionalGrid[j - 1])) * (S2 - Math.Cos(additionalGrid[j - 1]));

          b = 1.0 * ch_b / zn;
          c = 1.0 * ch_c / zn;
          a = 1 - b - c;
          return a * f(additionalGrid[j - 1]) + b * f(additionalGrid[j]) + c * f(additionalGrid[j + 1]);

        case ApproximationType.Hyperbolic:
          if (j == -2 || j == Configuration.GridPoints - 1)
          {
            return f(additionalGrid[j]);
          }

          S1 = 1.0 * (Math.Cosh(grid[j + 2]) - Math.Cosh(grid[j + 1])) / Math.Sinh(grid[j + 2] - grid[j + 1]);
          S2 = 1.0 * (Math.Sinh(grid[j + 2]) - Math.Sinh(grid[j + 1])) / Math.Sinh(grid[j + 2] - grid[j + 1]);

          zn = (Math.Sinh(additionalGrid[j + 1]) - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j]) - Math.Cosh(additionalGrid[j - 1]))
            - (Math.Sinh(additionalGrid[j]) - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j + 1]) - Math.Cosh(additionalGrid[j - 1]));

          ch_b = (Math.Sinh(additionalGrid[j + 1]) - Math.Sinh(additionalGrid[j - 1])) * (S2 - Math.Cosh(additionalGrid[j - 1]))
            - (S1 - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j + 1]) - Math.Cosh(additionalGrid[j - 1]));

          ch_c = (S1 - Math.Sinh(additionalGrid[j - 1])) * (Math.Cosh(additionalGrid[j]) - Math.Cosh(additionalGrid[j - 1]))
            - (Math.Sinh(additionalGrid[j]) - Math.Sinh(additionalGrid[j - 1])) * (S2 - Math.Cosh(additionalGrid[j - 1]));

          b = 1.0 * ch_b / zn;
          c = 1.0 * ch_c / zn;
          a = 1 - b - c;
          return a * f(additionalGrid[j - 1]) + b * f(additionalGrid[j]) + c * f(additionalGrid[j + 1]);

        default:
          throw new NotSupportedException($"Functional type {Configuration.ApproxType.ToString()} is not supported");
      }
    }
  }
}
