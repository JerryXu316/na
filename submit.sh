#!/bin/bash

# 设置你的 Git 仓库路径
REPO_PATH="$(pwd)"

# 添加所有更改到暂存区
git add .

# 提交更改
git commit -m "Submit homework"

# 推送到远程仓库的 main 分支
git push origin main

# 检查 Git 命令是否成功执行
if [ $? -eq 0 ]; then
    echo "Homework submitted successfully!"
else
    echo "Error: Failed to submit homework."
fi