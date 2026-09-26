# Statistical Moment Steering

MATLAB implementation of statistical moment steering for a nonlinear,
coupled-oscillator system. See <arxiv link>

## Requirements

The code was audited with MATLAB R2024b. The following MATLAB products are
used by the indicated workflows:

| Dependency | Used for |
| --- | --- |
| MATLAB | Dynamics propagation, optimization workflow, and plotting |
| Optimization Toolbox | `fsolve` in `generate_reference_curve.m` |
| Statistics and Machine Learning Toolbox | Chi-square quantiles, percentiles, fitted distributions, skewness, and kurtosis |
| Symbolic Math Toolbox | Variable-precision moment calculations |
| Parallel Computing Toolbox | The `parfor` loop in `main_MC.m` |
| [CVX with Mosek](https://cvxr.com/cvx/) | Convex subproblems in `main_SMS.m` and `main_CS.m` |

## Setup

Clone or download the repository, start MATLAB, and make the repository root
the current folder:

```matlab
cd("path/to/statistical-moment-steering")
```

Run all entry-point scripts from this directory. They call `initialization.m`
automatically, which adds `Supporting_Scripts` and its subfolders to the
MATLAB path.

## Workflow

The generated-data pipeline is:

```text
generate_reference_curve.m
        |
        v
     curve.mat
        |
        +-------------------+
        v                   v
    main_SMS.m          main_CS.m
        |                   |
        v                   v
 refTraj_*.mat        refTraj_CS.mat
        |                   |
        +---------+---------+
                  v
              main_MC.m
                  |
                  v
          results_MC_*.mat
```

### 1. Generate the reference curve

Run:

```matlab
generate_reference_curve
```

This uses a shooting method to determine the reference boundary conditions
and saves `curve.mat`. The remaining main scripts load this file through
`params_simulation.m`.

### 2. Generate an SMS solution with just mean and covariance constraint

Open `main_SMS.m` and select the objective:

```matlab
opt.objective = "quadratic";  % "quadratic", "cvar", or "dv99"
```

Disable the skewness constraint for the initial solve:

```matlab
opt.constraints.skew = 0;
```

Then run:

```matlab
main_SMS
```

The output is saved as `refTraj_<objective>.mat`, for example
`refTraj_quadratic.mat`.

### 3. Generate a skewness-constrained SMS solution

A base solution for the selected objective must already exist from Step 2.
Set:

```matlab
opt.constraints.skew = 1;
```

Run `main_SMS` again. It loads `refTraj_<objective>.mat` as the initial
reference and saves the result as `refTraj_<objective>_skew.mat`.

For the default quadratic configuration, the required sequence is therefore:

1. Run `main_SMS` with `opt.constraints.skew = 0` to create
   `refTraj_quadratic.mat`.
2. Run it again with `opt.constraints.skew = 1` to create
   `refTraj_quadratic_skew.mat`.

### 4. Generate the covariance-steering baseline

Run:

```matlab
main_CS
```

This solves the linearized covariance-steering problem and saves
`refTraj_CS.mat`.

### 5. Run Monte Carlo validation

Open `main_MC.m` and select a previously generated solution:

```matlab
solution_name = "quadratic";
```

Supported names correspond to available `refTraj_<solution_name>.mat` files,
such as `"quadratic"`, `"quadratic_skew"`, `"cvar"`, or `"CS"`.

Set the desired Monte Carlo indices. A small range is recommended for an
initial check:

```matlab
MC_indices = 0:100;
```

Then run:

```matlab
main_MC
```

The script calls `main_SIM.m` for each realization and saves
`results_MC_<solution_name>.mat`.

## Generated files

MATLAB data files are intentionally excluded by `.gitignore`. They must be
regenerated locally using the workflow above:

- `curve.mat`: reference boundary conditions and nonlinear model parameters
- `refTraj_*.mat`: optimized steering solutions and convergence histories
- `results_MC_*.mat`: Monte Carlo control-cost results
