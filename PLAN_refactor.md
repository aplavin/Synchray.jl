
### 2.2 Minimal refactor: make `synchrotron_model` explicit and extensible (no jet rewrite)

Motivation:

- Today, `AbstractSynchrotronMedium` is intentionally simple: geometry/flow are object methods, and microphysics enters only through
	- `electron_density(obj, x4)`
	- `magnetic_field(obj, x4)`
	- `synchrotron_model(obj)`
- Stage‑2 wants more variation (different spectra, different normalization semantics, later angle dependence), but we can get most of the value by *only* growing the `synchrotron_model` family first.

Goal:

- Allow swapping microphysics by changing `synchrotron_model(obj)` (and maybe its parameters), without touching jet geometry/kinematics code.
- Use an angle-capable transfer/microphysics entry point: `emissivity_absorption(obj, x4, k')`, where `k'` represents the photon in the plasma rest frame.

#### What stays the same

- No decomposition of jets into geometry/field/particles/process objects.
- No changes to `ConicalBKJet` / `ConicalBKJetWithPatterns` structure and delegation.
- The Stage‑1 invariant integration stays based on `(j_ν', α_ν')` computed in the comoving frame.
- Patterns remain microphysics-only modifiers (no geometry/flow changes).

#### Remaining work

Direction-aware microphysics plumbing (partial ordering + pitch-angle distributions)

True angle dependence needs (at least) the comoving photon direction $n'$ and the ordered field direction $b'$ (and later, electron pitch-angle distribution around $b'$). The key is to add this plumbing *without* refactoring jets into components.

**Current API shape**


- Medium entry point:
	- `emissivity_absorption(obj, x4, k')` is the primary/only interface.
	- `k'` is a proper photon 4-wavevector **in the plasma rest frame** (a null 4-vector).
		- The model extracts the comoving frequency $\nu'$ and direction $\hat n'$ directly from `k'`.
		- This avoids defining/transporting an arbitrary comoving spatial triad while still letting the model form invariants like $\mu=b'\cdot \hat n'$.

**Interpreting `k'` (what microphysics extracts)**

- Microphysics needs only (both are extracted from passed 4-vector `k'`):
	- comoving frequency $\nu'$
	- comoving propagation direction $\hat n'$ (unit 3-vector in the plasma rest frame)

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

- $\nu'$ and $\hat n'$ from `k'`
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

What the direction-aware `_synchrotron_coeffs` should receive (minimal):


- Extend the model dispatch for direction-aware models to accept the additional scalars needed:
	- `_synchrotron_coeffs(model, n_e, field, k')`
	- The model computes $\mu$ internally from the returned magnetic field and `k'`.

**Incremental rollout**

- First milestone: partially ordered/tangledness via $\kappa$.
- Second milestone: explicit pitch-angle distributions for electron anisotropy.

3) Keep “patterns” as they are

- `ConicalBKJetWithPatterns` remains the way to modulate `n_e` and `|B'|` without affecting geometry/flow.
