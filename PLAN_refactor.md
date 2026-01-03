
### 2.2 Minimal refactor: make `synchrotron_model` explicit and extensible (no jet rewrite)

Motivation:

- Today, `AbstractSynchrotronMedium` is intentionally simple: geometry/flow are object methods, and microphysics enters only through
	- `electron_density(obj, x4)`
	- `magnetic_field(obj, x4)`
	- `synchrotron_model(obj)`
- Stage‑2 wants more variation (different spectra, different normalization semantics, later angle dependence), but we can get most of the value by *only* growing the `synchrotron_model` family first.

Goal:

- Allow swapping microphysics by changing `synchrotron_model(obj)` (and maybe its parameters), without touching jet geometry/kinematics code.
- Move the transfer/microphysics entry point to an angle-capable representation: `emissivity_absorption(obj, x4, k')`, where `k'` represents the photon in the plasma rest frame.
	- This is a breaking change (replaces `emissivity_absorption(obj, x4, ν')`)

#### What stays the same

- No decomposition of jets into geometry/field/particles/process objects.
- No changes to `ConicalBKJet` / `ConicalBKJetWithPatterns` structure and delegation.
- The Stage‑1 invariant integration stays based on `(j_ν', α_ν')` computed in the comoving frame.

#### Minimal changes (recommended)

1) Rename the current Stage‑1 model to be explicit about its assumptions

- Rename `PowerLawElectrons` → `AngleAveragedPowerLawElectrons`.
	- Rationale: the current implementation *bakes in* the “isotropically tangled field” / angle‑averaged Stokes‑I assumption by averaging $\langle\sin^q\theta_{Bn}\rangle$ into the coefficients.

2) Introduce 1–2 new `synchrotron_model` types (still compatible with the current medium API)

These new models must work with the (updated) signature:

- `_synchrotron_coeffs(model, obj, x4, u, n_e, field, k')`

Suggested minimal set:

- `AngleAveragedPowerLawElectrons` (renamed current default; no behavior change).
- `FixedViewingAnglePowerLawElectrons(; p, θBn, ...)` (diagnostic model): treat $B_\perp = |B'|\sin\theta_{Bn}$ with user‑supplied $\theta_{Bn}$.
	- This model intentionally ignores the direction information in `k'` and is useful as a controlled sanity check.

3) Angle-dependent microphysics plumbing

True angle dependence needs (at least) the comoving photon direction $n'$ and the ordered field direction $b'$ (and later, electron pitch-angle distribution around $b'$). The key is to add this plumbing *without* refactoring jets into components.

**Minimal API shape**

- Replace the medium entry point with an angle-aware signature:
	- `emissivity_absorption(obj, x4, k')` is the primary/only interface.
	- `k'` is a proper photon 4-wavevector (a null 4-vector).
		- The comoving frequency is always computed as $\nu' = -k'\cdot u$.
		- The comoving propagation direction in the local rest space is derived as a rest-space vector
		  $$n' = k'/\nu' - u$$
		  with $u\cdot n' = 0$ and $n'\cdot n' = +1$.
		- This avoids defining/transporting an arbitrary comoving spatial triad while still letting the model form invariants like $\mu=b'\cdot n'$.

**Interpreting `k'` (what microphysics extracts)**

- Microphysics needs only (both are extracted from passed 4-vector `k'`):
	- comoving frequency $\nu'$
	- comoving propagation direction $n'$ (unit rest-space 4-vector; i.e. orthogonal to $u$)

**Ordered magnetic fields: what the microphysics must receive**


- For Stokes‑I with an *ordered* field direction, the radiative coefficients depend on
	- the comoving amplitude $|B'|$, and
	- the viewing geometry through $\mu \equiv b'\cdot n'$ (equivalently, $\theta_{Bn}$).
	- Typically $B_\perp = |B'|\,\sqrt{1-\mu^2}$.

Key design choice (to keep `_synchrotron_coeffs` signatures minimal): **$\mu$ is not passed as an argument**.

Instead, $\mu$ is computed inside the `synchrotron_model` from the returned magnetic field representation and the photon 4-vector `k'`.

Minimal representation choice (ties into the `magnetic_field(obj, x4)` rename):

- For the current Stage‑1 `AngleAveragedPowerLawElectrons`, have
	- `magnetic_field(obj, x4) -> FullyTangled(strength)`
	- i.e. return only the scalar comoving amplitude, wrapped to make the “tangled/angle‑averaged” assumption explicit.
- For new direction-aware models, allow `magnetic_field(obj, x4)` to return either:
	- `PartiallyTangled(field_vector; kappa)` where `field_vector` is the comoving ordered field vector (e.g. `SVector{3}`), or
	- just `field_vector::SVector` for fully ordered fields.

Then the model derives

- $\nu' = -k'\cdot u$
- $n' = k'/\nu' - u$ and $\hat n'$ from $n'$
- $\hat b'$ from the returned `field_vector`
- $\mu = \hat b'\cdot \hat n'$.

**Tangled/partially ordered fields**

- Introduce a minimal “ordering/tangledness” parameter $\kappa$.
	- For now, keep $\kappa$ constant per object/model (OK for performance and matches “same object → same return type”).
	- Natural convention (Fisher/concentration-style): $\kappa=0$ isotropic/tangled; $\kappa\to\infty$ fully ordered.
- The angle-dependent Stokes‑I model should consume only $(field, k', n_e)$, and compute $\nu'$, $\mu$, and the needed angle averages internally.
	- It can implement $\langle\sin^q\theta_{Bn}\rangle(\mu, \kappa)$ via lookup/caching inside `prepare_for_computations(model)` when that becomes performance-critical.

**Electron pitch-angle (velocity) distribution plumbing**

- Stage‑1 assumes isotropic electron momenta, so the pitch-angle distribution around $\hat b'$ is implicit.
- Angle-dependent and polarization-capable regimes eventually need an explicit model for electron pitch angles $\alpha_{eB}$:
	- Define a small family of pitch-angle distributions (configuration objects) as part of the `synchrotron_model`, e.g.
		- `IsotropicPitchAngles()` (status quo)
		- `FixedPitchAngle(α0)` (diagnostic)
		- `GaussianPitchAngle(α0, σα)`
	- The model then computes coefficients by integrating/averaging over $\alpha_{eB}$ as required by the chosen approximation.

What the angle-aware `_synchrotron_coeffs` should receive (minimal):

- Extend the model dispatch for angle-aware models to accept the additional scalars needed:
	- `_synchrotron_coeffs(model, obj, x4, u, n_e, field, k')`
	- The model computes $\mu$ internally from the returned magnetic field and `k'`.

**Incremental rollout**

- First milestone: ordered-field Stokes‑I with isotropic electrons.
	- Coefficients depend only on `(field, k', n_e)`; $\nu'$ and $\mu$ are computed internally.
- Second milestone: partially ordered/tangledness via $\kappa$.
- Third milestone: explicit pitch-angle distributions for electron anisotropy.

4) Keep “patterns” as they are

- `ConicalBKJetWithPatterns` remains the way to modulate `n_e` and `|B'|` without affecting geometry/flow.
