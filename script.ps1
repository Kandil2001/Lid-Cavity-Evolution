
git checkout main

# Temp folder outside the repo
$TempDir = "../temp_extract"
if (Test-Path $TempDir) {
    Remove-Item -Recurse -Force $TempDir
}
New-Item -ItemType Directory -Path $TempDir | Out-Null

# Create full directory structure
$folders = @(
    "matlab/iterative-solver",
    "matlab/vectorized-solver",
    "python/serial/iterative",
    "python/serial/vectorized"
)

foreach ($f in $folders) {
    New-Item -ItemType Directory -Path (Join-Path $TempDir $f) -Force | Out-Null
}

# --- MATLAB iterative-solver branch ---
git checkout iterative-solver
Copy-Item "IterativeSolver.m" (Join-Path $TempDir "matlab/iterative-solver") -Force
Copy-Item "README.md" (Join-Path $TempDir "matlab/iterative-solver") -Force

# --- MATLAB vectorized-solver branch ---
git checkout vectorized-solver
Copy-Item "VectorizedSolver.m" (Join-Path $TempDir "matlab/vectorized-solver") -Force
Copy-Item "README.md" (Join-Path $TempDir "matlab/vectorized-solver") -Force

# --- Python serial iterative branch ---
git checkout iterative
Copy-Item "IterativeSolver.py" (Join-Path $TempDir "python/serial/iterative") -Force
Copy-Item "README.md" (Join-Path $TempDir "python/serial/iterative") -Force

# --- Python serial vectorized branch ---
git checkout vectorized
Copy-Item "VectorizedSolver.py" (Join-Path $TempDir "python/serial/vectorized") -Force
Copy-Item "README.md" (Join-Path $TempDir "python/serial/vectorized") -Force

git checkout python

Copy-Item "README.md" (Join-Path $TempDir "python/") -Force
git checkout matlab

Copy-Item "README.md" (Join-Path $TempDir "matlab/") -Force



# Switch back to main
git checkout main

# Copy organized content into repo
Copy-Item (Join-Path $TempDir "*") . -Recurse -Force

# Stage & commit
git add .
git commit -m "new structure"

Write-Host "✅ Files copied into structured layout and committed on main."
