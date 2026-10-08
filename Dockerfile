FROM tercen/runtime-r44:4.4.3-11

WORKDIR /operator

# Tercen runs operators as UID 1000; default /root/.cache is unreadable.
ENV RENV_PATHS_CACHE=/operator/.cache/renv

# Dependencies first: this layer is only rebuilt when renv.lock changes,
# not on every code change.
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/settings.json renv/
# The base image already has the locked versions installed: copy them into the
# project library (hydrate) instead of compiling from source; restore() then
# only installs what differs from renv.lock.
RUN R -e "renv::consent(provided = TRUE);           renv::hydrate(packages = names(renv::lockfile_read()\$Packages),                         sources = c('/usr/local/lib/R/site-library', '/usr/local/lib/R/library'), prompt = FALSE);           renv::restore(confirm = FALSE)"

COPY . /operator
RUN chown -R 1000:1000 /operator

ENV TERCEN_SERVICE_URI=https://tercen.com

ENTRYPOINT ["R", "--no-save", "--no-restore", "--no-environ", "--slave", "-f", "main.R", "--args"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
