# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 07/11/2026
> __Description__: Create a laminar (plane Poiseuille) initial condition on the grid of sim_name:
>                      u = 1.5 U_b [1 - (y - y_mid)^2 / h^2],   v = w = 0     (outer units, U_b = 1)
>                  and save it in data/laminar_IC/ together with the simulation package files.
'''
# %%
# ~  LOAD THE FRAMEWORK
#~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
from io_operations import write_ucvcwc_file, read_ucvcwc_file
ts0 = tt.time()

# %%
# ~ POISEUILLE PROFILE ON THE STAGGERED GRID
# ~ ~~~~~~~~~~~~~~~~~~~~~
h     = 0.5 * (y_wall_top - y_wall_bot)                      #^ Channel half-height
y_mid = 0.5 * (y_wall_top + y_wall_bot)
U_yu  = 1.5 * (1.0 - ((yu - y_mid) / h)**2)                  #^ On u cell centers (with halos): ghost ~ -u_1 -> no-slip

u = np.broadcast_to(U_yu[None, None, :], (ngz + 1, ngx + 2, ngy + 1)).copy()   #^ [z, x, y]
v = np.zeros((ngz + 1, ngx + 1, ngy + 2))
w = np.zeros((ngz + 2, ngx + 1, ngy + 1))

dy_u  = np.diff(yy)                                          #^ Cell heights (interior u points: yu[1:-1])
U_b   = np.sum(U_yu[1:-1] * dy_u) / np.sum(dy_u)
log.info(f"Laminar profile: U_max = {U_yu.max():.4f}, U_b (discrete) = {U_b:.6f}")

# %%
# ~ SAVE
# ~ ~~~~~~~~~~~~~~~~~~~~~
save_loc = base_data + 'laminar_IC/'
create_directory(save_loc)

    # --> Copy the simulation package (grid, inputs, ...) except the velocity fields
vel_files = [pre + "{:07d}".format(frame_IC) for pre in ('01uc00.dble.', '02vc00.dble.', '03wc00.dble.')]
for fname in sorted(os.listdir(root_sim)):
    if fname not in vel_files and os.path.isfile(root_sim + fname):
        shutil.copy2(root_sim + fname, save_loc + fname)

    # --> Write and verify
write_ucvcwc_file(u, v, w, num=frame_IC, addr=save_loc)
_, u_chk, v_chk, w_chk = read_ucvcwc_file(save_loc, ngx, ngy, ngz, frame_IC)
assert np.array_equal(u_chk, u) and np.array_equal(v_chk, v) and np.array_equal(w_chk, w), "Read-back mismatch"

log.info(f"Laminar IC saved to: {save_loc}")
log.info(f"Total time taken is {tt.time() - ts0:.2f} seconds")
# %%
