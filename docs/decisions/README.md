# Architecture Decision Records

Each record describes one design decision: the context that forced it, the decision itself, and the consequences we accept, including the bad ones. Records are never rewritten after acceptance. A changed decision gets a new record that supersedes the old one.

| # | Decision | Status |
|---|---|---|
| [0001](0001-full-8d-phase-space.md) | Full 8D phase space and physics-agnostic integrators | Accepted |
| [0002](0002-numba-first.md) | Numba-first numerical core, SciPy as oracle, JAX optional | Accepted |
| [0003](0003-nfev-as-cost-metric.md) | Number of vector-field evaluations (nfev) as the primary cost metric | Accepted |
| [0004](0004-metric-interface.md) | Metric interface: compiled kernels returning full arrays | Accepted |

## Template

```markdown
# NNNN. Title

- Status: Proposed | Accepted | Superseded by NNNN
- Date: YYYY-MM-DD

## Context
What forces a decision, and what constraints apply.

## Decision
What we do, stated so that it can be checked in code review.

## Consequences
What becomes easier, what becomes harder, and what we must watch for.

## Alternatives considered
Each alternative, and why it was rejected.
```
