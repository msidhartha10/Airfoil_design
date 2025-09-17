
# NACA 4-Digit Airfoil Generator (MATLAB)

MATLAB implementation of the classical NACA 4-digit airfoil equations.  
Generates surface coordinates, plots the geometry, and exports results to Excel.  
Useful for CFD preprocessing, structural models, or visualization.

---

## ✈️ Features
- Generate any NACA 4-digit airfoil (e.g. `2412`, `5420`, `0015`)
- Cosine-spaced x-distribution (clusters points near LE/TE)
- Handles symmetric foils (`p = 0`)
- Outputs:
  - Upper/lower surface arrays
  - Camber line
  - Closed-loop surface array (TE → LE → TE)
- Built-in plotting
- Export to `.xlsx` for external tools

---

## 📖 Theory

The NACA 4-digit series defines an airfoil using four digits:

- 1st digit: maximum camber, as % of chord (`m`)
- 2nd digit: location of maximum camber from the LE, in tenths of chord (`p`)
- Last two digits: maximum thickness, as % of chord (`t`)

Example: **NACA 2412**
- Maximum camber = 2% of chord
- Camber location = 40% chord
- Thickness = 12% of chord

### Key Equations

**1. Thickness distribution**
```math
y_t = 5t \left( 0.2969 \sqrt{x} - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 + a_4 x^4 \right)
````

where
`a₄ = -0.1036` (closed TE), or `-0.1015` (open TE).

**2. Camber line**

```math
y_c =
\begin{cases}
\frac{m}{p^2}(2px - x^2), & x < p \\
\frac{m}{(1-p)^2}((1 - 2p) + 2px - x^2), & x ≥ p
\end{cases}
```

**3. Camber slope**

```math
\frac{dy_c}{dx} =
\begin{cases}
\frac{2m}{p^2}(p - x), & x < p \\
\frac{2m}{(1-p)^2}(p - x), & x ≥ p
\end{cases}
```

**4. Surface coordinates**

```math
\theta = \arctan\!\left(\frac{dy_c}{dx}\right)
```

```math
x_u = x - y_t \sin(\theta) \\
y_u = y_c + y_t \cos(\theta) \\
x_l = x + y_t \sin(\theta) \\
y_l = y_c - y_t \cos(\theta)
```

---

## 🚀 Usage

```matlab
c = 180;  % Chord length in mm

% Generate a NACA 5420 airfoil with 299 total points, closed trailing edge
[xsurf, ysurf, xu, yu, xl, yl, x, yc] = generate_naca4_airfoil('5420', c, 299, true);
```

This will:

1. Plot the NACA 5420 airfoil
2. Save coordinates as:

   ```
   NACA5420_airfoil_180mm.xlsx
   ```

---

## 📊 Example Plot

### NACA 5420 (Chord = 180 mm)

![Example Airfoil](docs/naca5420.png)

---

## 📂 File Structure

```
NACA_Airfoil_Generator/
├── generate_naca4_airfoil.m   % main function
├── example.m                  % demo usage script
├── README.md
└── LICENSE
```

---

## 📚 References

* Anderson, J. D. *Fundamentals of Aerodynamics*
* NASA Technical Reports on NACA airfoils
* Davies, A. J. — [Introduction to NACA Airfoil Aerodynamics in Python](https://towardsdatascience.com/introduction-to-naca-airfoil-aerodynamics-in-python-72a1c3ee46b1)

---

