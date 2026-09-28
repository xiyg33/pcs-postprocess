"""解析个人目录配置；所有相对路径均以配置文件所在目录为基准。"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


DEFAULT_PATHS_CONFIG = "postprocess_paths.local.json"


@dataclass(frozen=True)
class SourcePaths:
    results_root: Path
    input_subdir: str
    output_subdir: str
    raw_subdir: str

    def mode_dir(self, mode: str) -> Path:
        return self.results_root / f"results_{mode}"

    def input_dir(self, mode: str) -> Path:
        return self.mode_dir(mode) / self.input_subdir

    def output_dir(self, mode: str) -> Path:
        return self.mode_dir(mode) / self.output_subdir

    def raw_dir(self, mode: str) -> Path:
        return self.mode_dir(mode) / self.raw_subdir


@dataclass(frozen=True)
class PathsConfig:
    project_config: Path
    sources: dict[str, SourcePaths]

    @classmethod
    def load(cls, path: Path) -> "PathsConfig":
        path = Path(path).resolve()
        if not path.is_file():
            raise ValueError(f"Path config does not exist: {path}; copy postprocess_paths.example.json to {DEFAULT_PATHS_CONFIG}")
        raw = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"project_config", "sources"}:
            raise ValueError("Path config requires project_config and sources")
        if not isinstance(raw["project_config"], str) or not raw["project_config"].strip():
            raise ValueError("project_config must be a nonempty path")
        if not isinstance(raw["sources"], dict) or not raw["sources"]:
            raise ValueError("sources must be a nonempty object")

        def resolve(value: str) -> Path:
            candidate = Path(value).expanduser()
            return (candidate if candidate.is_absolute() else path.parent / candidate).resolve()

        sources = {}
        for name, entry in raw["sources"].items():
            if not isinstance(name, str) or not name.strip() or not isinstance(entry, dict):
                raise ValueError("Each source needs a name and an object of paths")
            if set(entry) - {"results_root", "input_subdir", "output_subdir", "raw_subdir"}:
                raise ValueError(f"Unknown path keys for source {name}")
            root = entry.get("results_root")
            if not isinstance(root, str) or not root.strip():
                raise ValueError(f"{name}.results_root must be a nonempty path")
            subdirs = {}
            for key, default in (("input_subdir", "canonical"),
                                 ("output_subdir", "postprocess"),
                                 ("raw_subdir", "hil_raw")):
                value = entry.get(key, default)
                # 子目录只允许单个目录名，避免写入另一模式或逃出结果根目录。
                if (not isinstance(value, str) or not value.strip() or
                        value in (".", "..") or "/" in value or "\\" in value):
                    raise ValueError(f"{name}.{key} must be a single directory name")
                subdirs[key] = value
            if subdirs["input_subdir"] == subdirs["output_subdir"]:
                raise ValueError(f"{name}: input_subdir and output_subdir must differ")
            sources[name] = SourcePaths(resolve(root), **subdirs)
        return cls(resolve(raw["project_config"]), sources)

    def source(self, name: str) -> SourcePaths:
        try:
            return self.sources[name]
        except KeyError as exc:
            raise ValueError(f"Unknown source {name!r}; configured sources: {', '.join(self.sources)}") from exc


def load_paths_config(path: Path | None = None) -> PathsConfig:
    if path is not None:
        return PathsConfig.load(path)
    cwd_config = Path.cwd() / DEFAULT_PATHS_CONFIG
    if cwd_config.is_file():
        return PathsConfig.load(cwd_config)
    # editable 安装可从任意当前目录找到仓库内的个人配置。
    package_config = Path(__file__).resolve().parents[2] / DEFAULT_PATHS_CONFIG
    return PathsConfig.load(package_config)
