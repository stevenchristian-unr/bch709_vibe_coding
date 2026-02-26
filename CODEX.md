# Codex Instructions for This Repo

## Purpose
This file defines how Codex should behave when working in this repository. Refer to it for all prompts and coding tasks in this folder.

## Core Behavior
- Prefer pragmatic, production-lean solutions over over-engineering.
- Keep changes minimal and scoped to the request.
- Ask a clarifying question only when a decision materially affects behavior or output.
- Default to safe, reversible operations.

## Repo Conventions
- Use `rg` for searching files and text.
- Do not use destructive commands (e.g., `rm -rf`, `git reset --hard`) unless explicitly asked.
- Respect existing code style and patterns.
- Keep comments short and only where they add clarity.

## Workflow Expectations
- Summarize what changed and where.
- List any tests you ran, or state if none were run.
- If changes are not possible or blocked, state why and propose the next step.

## Output Formatting
- Use concise Markdown.
- Avoid nested lists.
- Use inline code for paths and commands.

