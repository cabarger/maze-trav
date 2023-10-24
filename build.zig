const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const exe = b.addExecutable(.{
        .name = "maze-trav",
        .target = target,
        .optimize = optimize,
    });
    exe.addCSourceFiles(&[_][]const u8{
        "src/main.cpp",
        "src/cabarger_cs121_include.cpp",
    }, &[_][]const u8{});
    exe.linkSystemLibrary("ncurses");
    exe.linkLibCpp();

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);

    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);
}
