### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 56cba971-26e0-4835-9581-416e5b2bb8ba
# Loading required packages
using Rotations; using LinearAlgebra;

# ╔═╡ 7c26d6be-3afa-11f0-3af2-cb3c6834db75
md"""
# On the polar decomposition of the deformation gradient tensor

**Kevin Moerman**

The deformation gradient tensor $\mathbf{F}$ can map a line element in the initial state $\mathrm{d}\mathbf{X}$ to the corresponding line element (in terms of size and orientation) in the final state $\mathrm{d}\mathbf{x}$, which is expressable as: 

$$\mathrm{d}\mathbf{x}=\mathbf{F}\mathrm{d}\mathbf{X}$$

The polar decomposition of the deformation gradient tensor produces: 

$$\mathbf{F}=\mathbf{Q}\mathbf{U}=\mathbf{V}\mathbf{Q}$$

Where $\mathbf{Q}$ is a rotation tensor, and $\mathbf{U}$ and $\mathbf{V}$ are known as the the right and left stretch tensors respectively. Hence the deformation event due to $\mathbf{F}$ can be thought of as a stretching followed by a rotation ($\mathbf{Q}\mathbf{U}$) or equivalently a rotation followed by a stretching ($\mathbf{V}\mathbf{Q}$). Which in equation form can be presented as: 

$$\mathrm{d}\mathbf{x}=\rlap{\overbrace{\phantom{\mathbf{F}}}^{\mathbf{Q}\mathbf{U}}}\underbrace{\mathbf{F}}_{\mathbf{V}\mathbf{Q}}\mathrm{d}\mathbf{X}$$

Although the eigen vectors for $\mathbf{U}$ and $\mathbf{V}$ differ, they share the same eigen values, i.e. the principal stretches $\lambda_1$, $\lambda_2$, and $\lambda_3$.

Obtaining the rotation tensor $\mathbf{Q}$, and the stretch tensors $\mathbf{U}$ and $\mathbf{V}$, as well as their eigen vectors and eigen values, can be conveniently achieved through the use of the singular value decomposition of $\mathbf{F}$, which is given by: 

$$\mathbf{F}=\mathbf{W}\mathbf{\Sigma}\mathbf{R}$$

Here $\mathbf{W}$ and $\mathbf{R}$ are rotation tensors, and $\mathbf{\Sigma}$ is the singular value matrix. This enables the computation of the rotation tensor $\mathbf{Q}$ through: 

$$\mathbf{Q}=\mathbf{W}\mathbf{R}^\top$$

and the computation of the right stretch tensor $\mathbf{U}$ through: 

$$\mathbf{U}=\mathbf{R}^\top\mathbf{F}=\mathbf{R}\mathbf{\Sigma}\mathbf{R}^\top$$

and the computation of the left stretch tensor $\mathbf{V}$ through: 

$$\mathbf{V}=\mathbf{F}\mathbf{R}^\top=\mathbf{W}\mathbf{\Sigma}\mathbf{W}^\top$$

The eigenvalues or principal stretches are obtained from the singular value matrix $\mathbf{\Sigma}$ through: 

$$\lambda_i=\Sigma_{ii}$$

The eigenvectors for the right stretch tensor $\mathbf{U}$ (or more formally of $\mathbf{C}=\mathbf{F}^\top\mathbf{F}$ but with which they coincide) are obtained from the columns of $\mathbf{R}$: 

$$\mathbf{n_j}=\sum_{i=1}^3 R_{ij}\mathbf{e_i}$$

Similarly the eigenvectors for the left stretch tensor $\mathbf{V}$ (or more formally of $\mathbf{B}= \mathbf{F}^{} \mathbf{F}^\top$ but with which they coincide) are obtained from the columns of $\mathbf{W}$: 

$$\mathbf{m_j}=\sum_{i=1}^3 W_{ij}\mathbf{e_i}$$
"""

# ╔═╡ 84ce04a7-716f-46e5-948f-70d8fc6af880
md"""
## Numerical implementation
"""

# ╔═╡ ebad6d14-385a-4c9a-b583-d591dfa86063
md"""
Setting up required packages. We need some here for rotations and linear algebra. 
"""

# ╔═╡ 31a86696-ed77-4a79-b457-1e19562e9405
md"""
A simulated deformation gradient tensor $\mathbf{F}$ is now created. This is done by first defining some known principal stretches $\lambda_i$, using this to define a known right stretch tensor $\mathbf{U}$, and finally rotating this stretch tensor using a rotation matrix $\mathbf{Q}$ to obtain $\mathbf{F}$ from: 

$$\mathbf{F}=\mathbf{Q}\mathbf{U}$$
"""

# ╔═╡ 218327a9-94ec-4f6a-958a-a33d010aa34c
md"""
Let's first define the principal stretches $\lambda_i$: 
"""

# ╔═╡ 0691749b-d3b8-4787-a49a-b4f56faeed61
λ₁_true = 1.23

# ╔═╡ d840932f-f152-4fb7-9d52-6a12b09c9ca8
λ₂_true = 1.05

# ╔═╡ 44064371-9d73-4bae-b2f4-bbe14a44ad97
λ₃_true = 0.85

# ╔═╡ f1861957-a505-4dbb-b82b-209877fee5e3
md"""
Now use $\lambda_i$ to define $\mathbf{U}$, the right stretch tensor: 
"""

# ╔═╡ f3a03037-87d0-4404-b208-72f9a29f3bb1
U_true = Diagonal([λ₁_true, λ₂_true, λ₃_true])

# ╔═╡ a1e517c7-1db3-49ad-9515-122046a14252
md"""
Now let's create an arbitrary rotation tensor $\mathbf{Q}$, here a triplet of Euler angles, each $\frac{π}{4}$ in radians, is used: 
"""

# ╔═╡ 2d653777-5133-45e0-8fdb-c49135a3cfdf
# Create true rotation matrix Q
Q_true = RotXYZ(0.25*π,0.25*π,0.25*π)

# ╔═╡ 8502ba94-82d7-4013-bc25-411494559396
md"""
Finally this lets us define $\mathbf{F}$:
"""

# ╔═╡ c7bacbf9-8cf3-4e2e-93a0-94cd150271ae
# Create simulated deformation gradient tensor F
F = Q_true*U_true

# ╔═╡ a43de9d4-deb8-4d1d-b523-81200c59c480
md"""
Next our job is to retrieve $\mathbf{U}$, $\mathbf{V}$, $\mathbf{C}$, $\lambda_i$ (and various strain metrics) from $\mathbf{F}$. 
"""

# ╔═╡ 35e6ac2d-6da8-4aa7-bebb-edab3725bc43
md"""
## Using the singular value decomposition
One approach is to use the sinular value decomposition of $\mathbf{F}$, i.e.: 

$$\mathbf{F}=\mathbf{W}\mathbf{\Sigma}\mathbf{R}$$. 

In Julia this can be achieved using `LinearAlgebra`'s `svd` function. 
"""

# ╔═╡ 7c9759ee-e9d8-4b97-8eeb-14ccb1df7244


# ╔═╡ 2ebd78ba-7826-4466-b8f3-93a25e5b76d5
W, Σ, R = svd(F)

# ╔═╡ 17de4b76-f9f8-4ee1-8b0c-9bb526110730
Q_svd = W*R'

# ╔═╡ 9a1a36f4-b836-4e09-8b06-c5b3008e527f
U_svd = R'*F

# ╔═╡ b48e4062-a098-4d6c-becf-c8464b5ca724
V_svd = F*R'

# ╔═╡ a9c7bd65-a107-400a-a9a4-43a82bbfe399
λ_principal_svd = Σ

# ╔═╡ 15e55cbc-4d86-4b73-ba09-fc7c6f9b5b5d
n₁ = R*[1.0,0.0,0.0]

# ╔═╡ e73609c6-d0c8-4456-ba7d-1e1f16ea20b8
n₂ = R*[0.0,1.0,0.0]

# ╔═╡ 0a465eb1-4fba-4c41-bbf2-a0cd0f6277d4
n₃ = R*[0.0,0.0,1.0]

# ╔═╡ 685de779-7a3e-40cb-84c8-1fb5c60cf7e7
md"""
## Using the eigen decomposition of the right Cauchy-Green tensor
"""

# ╔═╡ d18399cb-e92a-4b61-9c33-6b6ab210ee3a
md"""
Alternatively the eigen decomposition can be used. First we obtain a symmetric tensor by computing the right Cauchy-Green tensor $\mathbf{C}$ from: 
$$\mathbf{C} = \mathbf{F}^{\top}\mathbf{F}$$
"""

# ╔═╡ a662c346-10cd-478e-a237-6d1bed1cdbd9
C = F'*F

# ╔═╡ d5179d55-8e1a-4bd1-a984-b2f7ec4acc54
md"""
The eigen values of $\mathbf{C}$ are the squared principle stretches ${\lambda_i}^2$
"""

# ╔═╡ 79d349e2-f51e-4f5f-b5c6-ce84a278bb79
λ_sq, Q_eig = eigen(C)

# ╔═╡ c4935cbc-3a12-42ed-a81b-fa18ed55e0a6
md"""
To obtain the principal stretches $\lambda_i$ we simply take the square root of the eigen values of $\mathbf{C}$
"""

# ╔═╡ 305b440c-575f-4c58-b2db-ce1c0db673ae
λ_principal = .√λ_sq

# ╔═╡ bf74da47-cae2-4267-a759-0ca6ef96b597
md"""
Next we can form the right stretch tensor $\mathbf{U}$ in the principal component coordinate system by forming a diagonal matrix using the principal stretches. 
"""

# ╔═╡ bdc9a490-5efa-417e-ae48-f520104821e8
U_prin_eig =  Diagonal(λ_principal)

# ╔═╡ 2462de8c-152b-4ae4-ae4b-a21f1a919b07
md"""
And the same can be done with the squared principal stretches to get $\mathbf{C}$ in the principal component system.
"""

# ╔═╡ 5e3aa452-b279-472f-a407-f288ce13ee71
C_prin_eig =  Diagonal(λ_sq)

# ╔═╡ e9e4db95-f1cc-4c67-9b04-bb95eb699a46
md"""
Next we can rotate the right stretch tensor in the principal component system to our regular system using the rotation tensor $\mathbf{Q}$
"""

# ╔═╡ fea755a3-6f62-4e69-9f8e-ee4d388fb972
U_eig = Q_eig*U_prin_eig*Q_eig'

# ╔═╡ d1882c78-163b-49ae-9d0b-e7e756acb089
md"""
And similarly for the right Cauchy green tensor. 
"""

# ╔═╡ 9ae10e35-a5bd-456e-9125-05a56adc1d91
C_eig = Q_eig*C_prin_eig*Q_eig'

# ╔═╡ ed5047ba-d9ef-4d43-8778-2a8a4d1c8193
md"""
## Deriving strain metrics
"""

# ╔═╡ 952020e1-f75c-4ea6-a876-760f7f7e5546
md"""
The natural or logarithmic strain: 

$$\mathbf{E}_{log}=log(\mathbf{U})$$
"""

# ╔═╡ 0c4f4b88-2621-4d8a-8762-dbb2ab4ad162
E_log_eig = Q_eig*Diagonal(log.(λ_principal))*Q_eig'

# ╔═╡ a56f0963-6cee-490a-a046-08f800296319
md"""
The Green-Lagrange strain: 

$$\mathbf{E}_{GL}=\frac{1}{2}(\mathbf{U}^2-\mathbf{I})$$
"""

# ╔═╡ 55674dfd-3c29-4d8f-8cab-81a7cda6518c
E_GL = Q_eig*Diagonal(0.5 .- (λ_principal.^2.0 .- 1.0))*Q_eig'

# ╔═╡ afb0ac48-8c5a-4f15-934f-d48c689152a4
md"""
The Hencky ("linear") strain: 

$$\mathbf{E}_{H}=\mathbf{U}-\mathbf{I}$$
"""

# ╔═╡ d74e5f03-7c99-4abd-a393-7b6479b15194
E_lin = Q_eig*Diagonal(λ_principal .- 1.0)*Q_eig'

# ╔═╡ ddec4711-b5e0-4871-a9c7-0f82681a2bbf
md"""
Some other Seth-Hill strain: 

$$\mathbf{E}_{SH}=\frac{1}{m}(\mathbf{U}^m-\mathbf{I})$$
"""

# ╔═╡ 451ef787-05c7-490e-a8f5-6437eccae5b2
m=3.0; E_SH = Q_eig*Diagonal(1.0/m .* (λ_principal.^m .- 1.0))*Q_eig'

# ╔═╡ a6f3a49b-683e-4373-ab16-0567aa95e28d


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Rotations = "6038ab10-8711-5258-84ad-4b1120ba62dc"

[compat]
Rotations = "~1.7.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "f8f02adba9601bc833e3f361b376696b51bf3dc4"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

    [deps.Rotations.weakdeps]
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─7c26d6be-3afa-11f0-3af2-cb3c6834db75
# ╟─84ce04a7-716f-46e5-948f-70d8fc6af880
# ╟─ebad6d14-385a-4c9a-b583-d591dfa86063
# ╠═56cba971-26e0-4835-9581-416e5b2bb8ba
# ╟─31a86696-ed77-4a79-b457-1e19562e9405
# ╟─218327a9-94ec-4f6a-958a-a33d010aa34c
# ╠═0691749b-d3b8-4787-a49a-b4f56faeed61
# ╠═d840932f-f152-4fb7-9d52-6a12b09c9ca8
# ╠═44064371-9d73-4bae-b2f4-bbe14a44ad97
# ╟─f1861957-a505-4dbb-b82b-209877fee5e3
# ╠═f3a03037-87d0-4404-b208-72f9a29f3bb1
# ╟─a1e517c7-1db3-49ad-9515-122046a14252
# ╠═2d653777-5133-45e0-8fdb-c49135a3cfdf
# ╟─8502ba94-82d7-4013-bc25-411494559396
# ╠═c7bacbf9-8cf3-4e2e-93a0-94cd150271ae
# ╟─a43de9d4-deb8-4d1d-b523-81200c59c480
# ╟─35e6ac2d-6da8-4aa7-bebb-edab3725bc43
# ╠═7c9759ee-e9d8-4b97-8eeb-14ccb1df7244
# ╠═2ebd78ba-7826-4466-b8f3-93a25e5b76d5
# ╠═17de4b76-f9f8-4ee1-8b0c-9bb526110730
# ╠═9a1a36f4-b836-4e09-8b06-c5b3008e527f
# ╠═b48e4062-a098-4d6c-becf-c8464b5ca724
# ╠═a9c7bd65-a107-400a-a9a4-43a82bbfe399
# ╠═15e55cbc-4d86-4b73-ba09-fc7c6f9b5b5d
# ╠═e73609c6-d0c8-4456-ba7d-1e1f16ea20b8
# ╠═0a465eb1-4fba-4c41-bbf2-a0cd0f6277d4
# ╟─685de779-7a3e-40cb-84c8-1fb5c60cf7e7
# ╟─d18399cb-e92a-4b61-9c33-6b6ab210ee3a
# ╠═a662c346-10cd-478e-a237-6d1bed1cdbd9
# ╠═d5179d55-8e1a-4bd1-a984-b2f7ec4acc54
# ╠═79d349e2-f51e-4f5f-b5c6-ce84a278bb79
# ╟─c4935cbc-3a12-42ed-a81b-fa18ed55e0a6
# ╠═305b440c-575f-4c58-b2db-ce1c0db673ae
# ╟─bf74da47-cae2-4267-a759-0ca6ef96b597
# ╠═bdc9a490-5efa-417e-ae48-f520104821e8
# ╟─2462de8c-152b-4ae4-ae4b-a21f1a919b07
# ╠═5e3aa452-b279-472f-a407-f288ce13ee71
# ╠═e9e4db95-f1cc-4c67-9b04-bb95eb699a46
# ╠═fea755a3-6f62-4e69-9f8e-ee4d388fb972
# ╟─d1882c78-163b-49ae-9d0b-e7e756acb089
# ╠═9ae10e35-a5bd-456e-9125-05a56adc1d91
# ╟─ed5047ba-d9ef-4d43-8778-2a8a4d1c8193
# ╟─952020e1-f75c-4ea6-a876-760f7f7e5546
# ╠═0c4f4b88-2621-4d8a-8762-dbb2ab4ad162
# ╟─a56f0963-6cee-490a-a046-08f800296319
# ╠═55674dfd-3c29-4d8f-8cab-81a7cda6518c
# ╟─afb0ac48-8c5a-4f15-934f-d48c689152a4
# ╠═d74e5f03-7c99-4abd-a393-7b6479b15194
# ╟─ddec4711-b5e0-4871-a9c7-0f82681a2bbf
# ╠═451ef787-05c7-490e-a8f5-6437eccae5b2
# ╠═a6f3a49b-683e-4373-ab16-0567aa95e28d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
