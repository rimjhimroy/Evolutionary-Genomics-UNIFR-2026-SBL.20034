# 1. Start from your base image
# Use the Renku VSCodium base image, which is compatible with Amalthea/Renku sessions.
FROM ghcr.io/swissdatasciencecenter/renku/py-basic-vscodium:2.15.0

ARG DEBIAN_FRONTEND=noninteractive

# 2. Switch to root to install packages
USER root

# Some base images may not ship with the expected non-root user.
# Create it if missing so later chown/USER/COPY --chown steps work.
ARG APP_USER=renku
RUN set -eux; \
        if id -u "$APP_USER" >/dev/null 2>&1; then \
            echo "$APP_USER already exists"; \
        else \
            if command -v useradd >/dev/null 2>&1; then \
                (command -v groupadd >/dev/null 2>&1 && groupadd -f "$APP_USER") || true; \
                useradd -m -s /bin/bash -g "$APP_USER" "$APP_USER"; \
            elif command -v adduser >/dev/null 2>&1; then \
                (command -v addgroup >/dev/null 2>&1 && addgroup -S "$APP_USER") || true; \
                adduser -S -G "$APP_USER" "$APP_USER"; \
            else \
                echo "Neither useradd nor adduser found; cannot create $APP_USER" >&2; \
                exit 1; \
            fi; \
        fi


# 3. Install dependencies (as root)
RUN mkdir -p /usr/share/man/man1 \
    && apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    curl \
    bzip2 \
    git \
    procps \
    openssl \
    locales \
    openjdk-17-jre-headless \
    fastqc \
    python3-pip \
    python3-venv \
    unzip \
    htop \
    && rm -rf /var/lib/apt/lists/* \
    # Install Quarto (for reproducible scientific documents)
    && curl -fsSL -o /tmp/quarto.deb https://quarto.org/download/latest/quarto-linux-amd64.deb \
    && apt-get update && apt-get install -y /tmp/quarto.deb \
    && rm /tmp/quarto.deb


# MultiQC: install into a dedicated venv (PEP 668 compliant)
RUN python3 -m venv /opt/multiqc-venv \
    && /opt/multiqc-venv/bin/pip install --no-cache-dir "multiqc==1.33" \
    && ln -sf /opt/multiqc-venv/bin/multiqc /usr/local/bin/multiqc

# 4. Optional: Add anything else you want to add as root
# Locale (helps R + some tools)
RUN sed -i 's/^# \(en_US.UTF-8 UTF-8\)/\1/' /etc/locale.gen \
    && locale-gen
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

COPY requirements.txt environment.yml /tmp/

## Pixi (pinned)
ENV PIXI_VERSION=0.65.0
ENV PIXI_TARBALL_SHA256=216adf73a7ba92547d38a4768c2c7cdad8cabb90af466b00c9ef69db89dc8649
RUN curl -fsSL -o /tmp/pixi.tar.gz \
      https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-x86_64-unknown-linux-musl.tar.gz \
    && echo "${PIXI_TARBALL_SHA256}  /tmp/pixi.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/pixi.tar.gz -C /usr/local/bin pixi \
    && chmod +x /usr/local/bin/pixi \
    && rm -f /tmp/pixi.tar.gz

# Install GATK + Picard (Java tools) directly
ENV GATK_VERSION=4.6.2.0
ENV PICARD_VERSION=3.4.0
RUN mkdir -p /opt/gatk /opt/picard \
    && curl -fsSL -o /tmp/gatk.zip \
      https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
    && unzip -q /tmp/gatk.zip -d /opt/gatk \
    && rm -f /tmp/gatk.zip \
    && ln -sf /opt/gatk/gatk-${GATK_VERSION}/gatk /usr/local/bin/gatk \
    && curl -fsSL -o /opt/picard/picard.jar \
      https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
    && printf '%s\n' '#!/usr/bin/env bash' 'exec java -jar /opt/picard/picard.jar "$@"' > /usr/local/bin/picard \
    && chmod +x /usr/local/bin/picard

RUN mkdir -p /opt/practical1 /work \
    && chown -R renku:renku /opt/practical1 /work


# 5. Switch back to the original user
USER renku

# 6. Optional: Add anything else you want to add as user
# Install Practical 1 environment (from lockfile)
WORKDIR /opt/practical1
COPY --chown=renku:renku pixi.toml pixi.lock /opt/practical1/
RUN pixi install --locked --no-progress

# Make tools available without "pixi run"
ENV PATH=/opt/practical1/.pixi/envs/default/bin:$PATH

# Seed VS Code defaults for Renku sessions.
COPY --chown=renku:renku .vscode /opt/practical1/.vscode
COPY --chown=renku:renku vscode-user /opt/practical1/vscode-user
RUN set -eux; \
    for user_dir in \
        /home/renku/.config/VSCodium/User \
        "/home/renku/.config/Code - OSS/User"; do \
        mkdir -p "$user_dir"; \
        cp /opt/practical1/vscode-user/settings.json "$user_dir/settings.json"; \
        cp /opt/practical1/vscode-user/keybindings.json "$user_dir/keybindings.json"; \
    done

# Include the tutorial README in the image
COPY --chown=renku:renku README_Practical1.md /opt/practical1/README_Practical1.md
COPY --chown=renku:renku README_Practical1.md /home/renku/work/README_Practical1.md


# IMPORTANT: do not override the base image CMD/ENTRYPOINT.
# In Kubernetes (Renku/Amalthea sessions) a non-interactive `bash` CMD will exit immediately,
# causing CrashLoopBackOff.
