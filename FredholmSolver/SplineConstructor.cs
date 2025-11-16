using System;
using System.Collections.Generic;

namespace FredholmSolver
{
  public class SplineConstructor
  {
    /// <summary>
    /// Построение значения минимального сплайна в точке.
    /// Значение определяется типом и номером сплайн-функции, а также расчётной сеткой.
    /// </summary>
    public double SplineValueAtPoint(   
      int splineNumber,
      Dictionary<int, double> grid,
      double point     
      )
    {
      switch (Configuration.ApproxType)
      {
        case ApproximationType.ShoenbergMarsden:
        case ApproximationType.Averaging:
        case ApproximationType.DeBoorFix0:
        case ApproximationType.DeBoorFix1:
        case ApproximationType.DeBoorFix2:
        case  ApproximationType.Proectional:
          return BSplineAtPoint(grid, splineNumber, point);
        case ApproximationType.AveragingTrigonometric:
        case ApproximationType.ProectionalTrigonometric:
        case ApproximationType.DeBoorFix1Trigonometric:
        case ApproximationType.DeBoorFix2Trigonometric:
          return TrigSplineAtPoint(grid, splineNumber, point);
        case ApproximationType.AveragingHyperbolic:
        case ApproximationType.ProectionalHyperbolic:
          return HypSplineAtPoint(grid, splineNumber, point);
        default:
          throw new NotSupportedException($"Spline type {Configuration.ApproxType.ToString()} is not supported");
      }
    }

    private double BSplineAtPoint(Dictionary<int, double> grid, int splineNumber, double point)
    {
      int j = splineNumber;
      if (point >= grid[j] && point < grid[j + 1])
      {
        double ch = (point - grid[j]) * (point - grid[j]);
        double zn = (grid[j] - grid[j + 1]) * (grid[j] - grid[j + 2]);

        return 1.0 * ch / zn;
      }
      else if (point >= grid[j + 1] && point < grid[j + 2])
      {
        double coeff = 1.0 / (grid[j + 1] - grid[j]);
        double first = 1.0 * (point - grid[j]) * (point - grid[j]) / (grid[j + 2] - grid[j]);
        double secondCh = (point - grid[j + 1]) * (point - grid[j + 1]) * (grid[j + 3] - grid[j]);
        double secondZn = (grid[j + 2] - grid[j + 1]) * (grid[j + 3] - grid[j + 1]);
        double second = 1.0 * secondCh / secondZn;

        return coeff * (first - second);
      }
      else if (point >= grid[j + 2] && point < grid[j + 3])
      {
        double ch = (point - grid[j + 3]) * (point - grid[j + 3]);
        double zn = (grid[j + 3] - grid[j + 1]) * (grid[j + 3] - grid[j + 2]);

        return 1.0 * ch / zn;
      }
      else
      {
        return 0;
      }
    }

    private double TrigSplineAtPoint(Dictionary<int, double> grid, int splineNumber, double point)
    {
      int j = splineNumber;
      if (point >= grid[j] && point < grid[j + 1])
      {
        double ch = Math.Cos(0.5 * (grid[j + 2] - grid[j + 1])) * 
          Math.Sin(0.5 * (point - grid[j])) * 
          Math.Sin(0.5 * (point - grid[j]));

        double zn = Math.Sin(0.5 * (grid[j + 1] - grid[j])) * 
          Math.Sin(0.5 * (grid[j + 2] - grid[j]));

        return 1.0 * ch / zn;
      }
      else if (point >= grid[j + 1] && point < grid[j + 2])
      {
        double coeff = 1.0 * Math.Cos(0.5 * (grid[j + 2] - grid[j + 1])) / 
          Math.Sin(0.5 * (grid[j + 1] - grid[j]));
        double first = 1.0 * Math.Sin(0.5 * (point - grid[j])) * 
          Math.Sin(0.5 * (point - grid[j])) / 
          Math.Sin(0.5 * (grid[j + 2] - grid[j]));

        double secondCh = Math.Sin(0.5 * (point - grid[j + 1])) * 
          Math.Sin(0.5 * (point - grid[j + 1])) * 
          Math.Sin(0.5 * (grid[j + 3] - grid[j]));
        double secondZn = Math.Sin(0.5 * (grid[j + 3] - grid[j + 1])) * 
          Math.Sin(0.5 * (grid[j + 2] - grid[j + 1]));
        double second = 1.0 * secondCh / secondZn;

        return coeff * (first - second);
      }
      else if (point >= grid[j + 2] && point < grid[j + 3])
      {
        double ch = Math.Cos(0.5 * (grid[j + 2] - grid[j + 1])) * 
          Math.Sin(0.5 * (grid[j + 3] - point)) *
          Math.Sin(0.5 * (grid[j + 3] - point));
        double zn = Math.Sin(0.5 * (grid[j + 3] - grid[j + 1])) * 
          Math.Sin(0.5 * (grid[j + 3] - grid[j + 2]));

        return 1.0 * ch / zn;
      }
      else
      {
        return 0;
      }
    }

    private double HypSplineAtPoint(Dictionary<int, double> grid, int splineNumber, double point)
    {
      int j = splineNumber;
      if (point >= grid[j] && point < grid[j + 1])
      {
        double ch = Math.Cosh(0.5 * (grid[j + 2] - grid[j + 1])) * 
          Math.Sinh(0.5 * (point - grid[j])) * 
          Math.Sinh(0.5 * (point - grid[j]));
        double zn = Math.Sinh(0.5 * (grid[j + 1] - grid[j])) * 
          Math.Sinh(0.5 * (grid[j + 2] - grid[j]));

        return 1.0 * ch / zn;
      }
      else if (point >= grid[j + 1] && point < grid[j + 2])
      {
        double coeff = 1.0 * Math.Cosh(0.5 * (grid[j + 2] - grid[j + 1])) / 
          Math.Sinh(0.5 * (grid[j + 1] - grid[j]));
        double first = 1.0 * Math.Sinh(0.5 * (point - grid[j])) * 
          Math.Sinh(0.5 * (point - grid[j])) / 
          Math.Sinh(0.5 * (grid[j + 2] - grid[j]));

        double secondCh = Math.Sinh(0.5 * (point - grid[j + 1])) * 
          Math.Sinh(0.5 * (point - grid[j + 1])) * 
          Math.Sinh(0.5 * (grid[j + 3] - grid[j]));
        double secondZn = Math.Sinh(0.5 * (grid[j + 3] - grid[j + 1])) * 
          Math.Sinh(0.5 * (grid[j + 2] - grid[j + 1]));
        double second = 1.0 * secondCh / secondZn;

        return coeff * (first - second);
      }
      else if (point >= grid[j + 2] && point < grid[j + 3])
      {
        double ch = Math.Cosh(0.5 * (grid[j + 2] - grid[j + 1])) * 
          Math.Sinh(0.5 * (grid[j + 3] - point)) * 
          Math.Sinh(0.5 * (grid[j + 3] - point));
        double zn = Math.Sinh(0.5 * (grid[j + 3] - grid[j + 1])) * 
          Math.Sinh(0.5 * (grid[j + 3] - grid[j + 2]));

        return 1.0 * ch / zn;
      }
      else
      {
        return 0;
      }
    }
  }
}
