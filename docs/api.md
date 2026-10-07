# API reference

The package is `geoint`. Its layers depend on each other in one direction only:
`metrics → formulations → integrators (no physics) → testcases → experiments → plotting`.

## Metrics

::: geoint.metrics.base

::: geoint.metrics.schwarzschild.Schwarzschild

## Formulations

::: geoint.formulations.base.Formulation

::: geoint.formulations.second_order.SecondOrder

::: geoint.formulations.hamiltonian.Hamiltonian

## Initial data

::: geoint.initial_conditions

## Integrators

::: geoint.integrators.base

::: geoint.integrators.registry

## Analytic references

::: geoint.analytic.schwarzschild_deflection

::: geoint.analytic.schwarzschild_orbits

## Test cases

::: geoint.testcases.base

::: geoint.testcases.schwarzschild

## Experiments

::: geoint.experiments.runner

::: geoint.experiments.longterm

::: geoint.experiments.analysis
