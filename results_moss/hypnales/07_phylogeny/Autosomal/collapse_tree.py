import glob
from ete3 import Tree

# 1. 현재 디렉터리의 모든 .treefile 목록 가져오기
tree_files = glob.glob("G*.treefile")

# 임계값 설정 (필요에 따라 10, 33, 50, 70 등으로 수정)
THRESHOLD = 10
print(f"총 {len(tree_files)}개의 트리를 처리합니다...\n")

for tree_file in tree_files:
    t = Tree(tree_file)
    
    # 2. 트리의 모든 노드를 순회하며 지지도 확인
    for node in t.get_descendants():
        # 잎(종 이름)이 아니고, 부트스트랩 점수가 THRESHOLD 미만인 경우
        if not node.is_leaf() and node.support < THRESHOLD:
            # 해당 노드를 삭제하고 자식들을 부모에게 끌어올려 다분기(Polytomy) 형성
            node.delete()
            
    # 3. 처리된 트리를 새로운 파일명으로 저장
    output_name = "10_collapsed_" + tree_file
    t.write(outfile=output_name)
    print(f"처리 완료: {tree_file} -> {output_name}")
