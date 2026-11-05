import re, sys
src = open('verify.cpp').read()
def rep(old, new, count=1):
    global src
    assert src.count(old) == count, (old[:60], src.count(old))
    src = src.replace(old, new)
rep('enum class ExitCode : std::uint8_t { Success = 0U, Unresolved = 1U, UsageError = 2U, InputError = 3U, LimitError = 5U };',
    'enum class ExitCode : std::uint8_t { Success = 0U, Unresolved = 1U, UsageError = 2U, InputError = 3U, LimitError = 5U, Checkpoint = 7U };')
rep('''    Unresolved,
    StackExhausted
};''', '''    Unresolved,
    StackExhausted,
    Checkpoint
};''')
rep('''    case Error::StackExhausted:
        return "the stack of boxes is full";''', '''    case Error::StackExhausted:
        return "the stack of boxes is full";
    case Error::Checkpoint:
        return "node budget reached, checkpoint written";''')
rep('''    [[nodiscard]] std::expected<void, Error> run() {
        const Parameters& parameters{certificate_.parameters};
        const Value scale{parameters.scale};
        const Wide spanTarget{Wide{certificate_.rewardSumDouble + certificate_.potentialBudget} * scale};
        statistics_ = Statistics{};
        stack_[0U].lower.fill(parameters.minimumGap);
        stack_[0U].upper.fill(parameters.gridLength);
        std::size_t stackSize{1U};
        statistics_.maximumStack = 1U;
        while (stackSize > 0U) {
            --stackSize;''', '''    // N8OPT PATCH: checkpoint/resume. With an empty resume list the search starts from the full box
    // (as the reference program); otherwise it continues from the saved boxes and statistics. After
    // budget boxes it stops and keeps the waiting boxes (savedStack) for the next run.
    [[nodiscard]] std::expected<void, Error> run(const std::uint64_t budget, const std::vector<Box>& resume,
                                                 const Statistics& resumeStatistics) {
        const Parameters& parameters{certificate_.parameters};
        const Value scale{parameters.scale};
        const Wide spanTarget{Wide{certificate_.rewardSumDouble + certificate_.potentialBudget} * scale};
        std::size_t stackSize{0U};
        if (resume.empty()) {
            statistics_ = Statistics{};
            stack_[0U].lower.fill(parameters.minimumGap);
            stack_[0U].upper.fill(parameters.gridLength);
            stackSize = 1U;
            statistics_.maximumStack = 1U;
        } else {
            if ((resume.size() + 2U) >= kStackCapacity) {
                return std::unexpected(Error::StackExhausted);
            }
            statistics_ = resumeStatistics;
            for (const Box& saved : resume) {
                stack_[stackSize] = saved;
                ++stackSize;
            }
        }
        std::uint64_t processed{0U};
        while (stackSize > 0U) {
            if (processed >= budget) {
                savedStackSize_ = stackSize;
                return std::unexpected(Error::Checkpoint);
            }
            ++processed;
            --stackSize;''')
rep('''    [[nodiscard]] const Box& lastBox() const noexcept {
        return lastBox_;
    }''', '''    [[nodiscard]] const Box& lastBox() const noexcept {
        return lastBox_;
    }

    [[nodiscard]] std::span<const Box> savedStack() const noexcept {
        return std::span<const Box>{stack_}.first(savedStackSize_);
    }''')
rep('''    std::array<Box, kStackCapacity> stack_{};
};''', '''    std::array<Box, kStackCapacity> stack_{};
    std::size_t savedStackSize_{0U};
};''')
# fingerprint of the parameters (all integers in order) to tie checkpoints to one certificate
rep('''[[nodiscard]] ExitCode run(const std::span<const std::string> arguments) {
    if (arguments.size() != 3U) {
        std::cerr << "usage: verify params.txt kernel_lower.bin\\n";
        return ExitCode::UsageError;
    }''', '''[[nodiscard]] std::uint64_t fingerprint(const Parameters& p) noexcept {
    constexpr std::uint64_t kMod{(std::uint64_t{1} << 61U) - 1U};
    std::uint64_t h{12345U};
    const auto mix = [&h](const Value v) noexcept {
        const std::uint64_t u{static_cast<std::uint64_t>(v + (Value{1} << 40U))};
        __extension__ using U128 = unsigned __int128;
        h = static_cast<std::uint64_t>(((U128{h} * 1000003U) + u) % kMod);
    };
    mix(p.gridLength); mix(p.gridPerUnit); mix(p.scale); mix(p.alpha); mix(p.minimumGap); mix(p.firstKnot); mix(p.knotStep);
    for (const auto& row : p.weights) { for (const Value v : row) { mix(v); } }
    for (const auto& row : p.rewards) { for (const Value v : row) { mix(v); } }
    for (const auto& gap : p.potentials) { for (const auto& row : gap) { for (const Value v : row) { mix(v); } } }
    return h;
}

[[nodiscard]] ExitCode run(const std::span<const std::string> arguments) {
    if ((arguments.size() != 3U) && (arguments.size() != 5U) && (arguments.size() != 6U)) {
        std::cerr << "usage: verify params.txt kernel_lower.bin [budget checkpoint_out [checkpoint_in]]\\n";
        return ExitCode::UsageError;
    }''')
rep('''    BoxSearch search{certificate.value(), std::move(table).value()};
    const std::expected<void, Error> outcome{search.run()};
    if (!outcome.has_value()) {''', '''    const std::uint64_t print{fingerprint(certificate.value().parameters)};
    std::uint64_t budget{std::numeric_limits<std::uint64_t>::max()};
    std::vector<Box> resume{};
    Statistics resumeStatistics{};
    if (arguments.size() >= 5U) {
        budget = std::stoull(arguments[3U]);
    }
    if (arguments.size() == 6U) {
        std::ifstream in{arguments[5U]};
        std::string tag{};
        std::uint64_t version{0U};
        std::uint64_t savedPrint{0U};
        std::size_t count{0U};
        in >> tag >> version >> savedPrint >> resumeStatistics.nodes >> resumeStatistics.leaves >> resumeStatistics.basic >>
            resumeStatistics.affine >> resumeStatistics.maximumStack >> count;
        if (!in || (tag != "checkpoint") || (version != 1U) || (savedPrint != print) || (count == 0U) ||
            ((count + 3U) >= kStackCapacity)) {
            std::cerr << "verify: bad checkpoint\\n";
            return ExitCode::InputError;
        }
        const Parameters& p{certificate.value().parameters};
        for (std::size_t k = 0U; k < count; ++k) {
            Box b{};
            for (std::size_t gap = 0U; gap < kGapCount; ++gap) {
                in >> b.lower[gap] >> b.upper[gap];
                if (!in || (b.lower[gap] < p.minimumGap) || (b.upper[gap] > p.gridLength) || (b.lower[gap] > b.upper[gap])) {
                    std::cerr << "verify: bad box in checkpoint\\n";
                    return ExitCode::InputError;
                }
            }
            resume.push_back(b);
        }
        std::string endTag{};
        in >> endTag;
        if (endTag != "end") {
            std::cerr << "verify: checkpoint not terminated\\n";
            return ExitCode::InputError;
        }
    }
    BoxSearch search{certificate.value(), std::move(table).value()};
    const std::expected<void, Error> outcome{search.run(budget, resume, resumeStatistics)};
    if (!outcome.has_value() && (outcome.error() == Error::Checkpoint)) {
        std::ofstream out{arguments[4U]};
        const Statistics& s{search.statistics()};
        out << "checkpoint 1 " << print << ' ' << s.nodes << ' ' << s.leaves << ' ' << s.basic << ' ' << s.affine << ' '
            << s.maximumStack << ' ' << search.savedStack().size() << '\\n';
        for (const Box& b : search.savedStack()) {
            for (std::size_t gap = 0U; gap < kGapCount; ++gap) {
                out << b.lower[gap] << ' ' << b.upper[gap] << ' ';
            }
            out << '\\n';
        }
        out << "end\\n";
        out.close();
        if (!out) {
            std::cerr << "verify: cannot write checkpoint\\n";
            return ExitCode::LimitError;
        }
        writeStatistics(std::cout, "CHECKPOINT", s);
        return ExitCode::Checkpoint;
    }
    if (!outcome.has_value()) {''')
open('verify_cp.cpp', 'w').write(src)
print('patched ok')
